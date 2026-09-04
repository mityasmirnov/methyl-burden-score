#!/usr/bin/env bash
# Matched 16-epoch promotion screen (7G′ Stage A).
# N-light-gene-max 16-ep is COMPLETE (3/3) — do not rerun. Remaining = 12
# GPU-fold jobs (light-mean + scalar×2 + vector mean-max) on GPU 0.
# Prefer scripts/run_7g_16ep_promotion_resume.sh while the queue is in flight.
# See docs/plans/milestone-7g-prime-16ep-promotion.md
set -euo pipefail
cd "$(dirname "$0")/.."
source scripts/activate_data_environment.sh

export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
export PYTHONUNBUFFERED=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

MAIN_CFG=configs/experiment/stage0_7g_gene_only_probe.yaml
RUN=scripts/run_7g_gene_only_probe.py
REPORT=scripts/write_7g_gene_only_probe_report.py
ART=artifacts/runs
LOG_DIR=scratch/logs
mkdir -p "$LOG_DIR"
LOG="$LOG_DIR/7g_16ep_promotion_$(date +%Y%m%d_%H%M%S).log"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

# Abort if GPU 0 (visible as cuda:0) already has a compute process using >2 GiB.
require_gpu_free() {
  local used
  used=$(nvidia-smi --id=0 --query-gpu=memory.used --format=csv,noheader,nounits 2>/dev/null | head -1 | tr -d ' ')
  if [[ -z "$used" ]]; then
    log "ERROR: nvidia-smi failed for GPU 0"
    exit 1
  fi
  if (( used > 4096 )); then
    log "ERROR: GPU 0 has ${used} MiB in use — aborting before next fold"
    nvidia-smi --id=0 || true
    exit 1
  fi
  log "GPU 0 free enough (${used} MiB used)"
}

sync_flat_arm() {
  local arm_id="$1"
  local prefix="$2"
  python3 - <<PY
import json
from pathlib import Path
art = Path("$ART")
report = Path("reports/inspection/stage0_7g_gene_only_probe/per_arm")
report.mkdir(parents=True, exist_ok=True)
folds = []
for i in range(3):
    p = art / f"$prefix-f{i}" / "metrics.json"
    if p.is_file():
        folds.append(json.loads(p.read_text()))
out = {"arm_id": "$arm_id", "kind": "flat_region_train", "folds": folds, "n_folds": len(folds)}
(report / f"$arm_id.json").write_text(json.dumps(out, indent=2, default=str), encoding="utf-8")
print(f"synced $arm_id n_folds={len(folds)}")
PY
  uv run python "$REPORT" 2>&1 | tee -a "$LOG"
}

sync_cascade_arm() {
  local arm_id="$1"
  local run_id="$2"
  python3 - <<PY
import json
from pathlib import Path
art = Path("$ART") / "$run_id"
report = Path("reports/inspection/stage0_7g_gene_only_probe/per_arm")
report.mkdir(parents=True, exist_ok=True)
folds = []
for fold_dir in sorted(art.glob("fold_*")):
    p = fold_dir / "metrics.json"
    if p.is_file():
        folds.append(json.loads(p.read_text()))
out = {"arm_id": "$arm_id", "kind": "cascade_train", "run_id": "$run_id", "folds": folds, "n_folds": len(folds)}
(report / f"$arm_id.json").write_text(json.dumps(out, indent=2, default=str), encoding="utf-8")
print(f"synced $arm_id n_folds={len(folds)}")
PY
  uv run python "$REPORT" 2>&1 | tee -a "$LOG"
}

posthoc_nested_cascade() {
  local run_id="$1"
  log "SKIP nested enet for $run_id (post-hoc later: eval_mbs_enet_from_scores.py --nested)"
}

posthoc_nested_flat() {
  local prefix="$1"
  log "SKIP nested enet for $prefix (post-hoc later: eval_mbs_enet_from_scores.py --nested)"
}

log "=== 16-epoch promotion screen start (CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES) ==="
require_gpu_free

# N-light-gene-max 16-ep complete (3/3) — sync only, do not retrain.
log "=== N-light-gene-max already complete (skip train) ==="
sync_flat_arm "N-light-gene-max" "stage0-7g-gene-probe-light-max"
posthoc_nested_flat "stage0-7g-gene-probe-light-max"

# Remaining 12 folds: mean + scalar×2 + vector mean-max.
log "=== N-light-gene-mean all folds ==="
require_gpu_free
uv run python "$RUN" --config "$MAIN_CFG" --device cuda --arm "N-light-gene-mean" 2>&1 | tee -a "$LOG"
sync_flat_arm "N-light-gene-mean" "stage0-7g-gene-probe-light-mean-16ep"
posthoc_nested_flat "stage0-7g-gene-probe-light-mean-16ep"

log "=== N-cascade-scalar-mean-max ==="
require_gpu_free
uv run python "$RUN" --config "$MAIN_CFG" --device cuda --arm "N-cascade-scalar-mean-max" 2>&1 | tee -a "$LOG"
posthoc_nested_cascade "stage0-7g-gene-probe-scalar-mean-max-16ep"
sync_cascade_arm "N-cascade-scalar-mean-max" "stage0-7g-gene-probe-scalar-mean-max-16ep"

log "=== N-cascade-scalar-max-mean ==="
require_gpu_free
uv run python "$RUN" --config "$MAIN_CFG" --device cuda --arm "N-cascade-scalar-max-mean" 2>&1 | tee -a "$LOG"
posthoc_nested_cascade "stage0-7g-gene-probe-scalar-max-mean-16ep"
sync_cascade_arm "N-cascade-scalar-max-mean" "stage0-7g-gene-probe-scalar-max-mean-16ep"

log "=== N-cascade-vector-mean-max ==="
require_gpu_free
uv run python "$RUN" --config "$MAIN_CFG" --device cuda --arm "N-cascade-vector-mean-max" 2>&1 | tee -a "$LOG"
posthoc_nested_cascade "stage0-7g-gene-probe-vector-mean-max-16ep"
sync_cascade_arm "N-cascade-vector-mean-max" "stage0-7g-gene-probe-vector-mean-max-16ep"

log "=== nested enet SKIPPED (run post-hoc after GPU queue) ==="
log "  See scripts/run_7g_16ep_promotion_resume.sh footer for exact eval_mbs_enet_from_scores.py commands."

log "=== final report refresh ==="
uv run python "$REPORT" 2>&1 | tee -a "$LOG"
uv run python scripts/apply_7g_16ep_decision.py 2>&1 | tee -a "$LOG"
log "=== 16-epoch promotion screen complete ==="
log "Apply decision rules in reports/inspection/stage0_7g_gene_only_probe/promotion_decision.json"
