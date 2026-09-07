#!/usr/bin/env bash
# Milestone 10 overnight GPU-0 chain (2026-09-07 -> 2026-09-08 09:30 CEST).
#
# Replaces the first version of this script, which had a self-matching
# pgrep bug: `pgrep -f 'run_7h_next_queue.sh|...'` matches the *text* of
# any command line containing that filename, including unrelated diagnostic
# commands that merely mention it -- it doesn't require the process to be
# a real instance of that script. That kept this chain stuck "waiting"
# for over an hour after the thing it was waiting for had already crashed
# and exited, wasting GPU-0 time the user explicitly asked to keep busy.
# This version never text-matches process lists; it runs each step
# directly, in the foreground, one after another -- correctness by
# construction rather than by pattern-matching other processes.
#
# scripts/run_7h_next_queue.sh crashed at 19:57 CEST on an import bug in
# run_7h_onehop_correctness_smokes.py (open_betas_for_matrix lives in
# mbs.matrix.virtual_hub_store, not mbs.matrix.store) -- fixed separately.
# That crash meant Track B.4 (queued after the smokes) never ran either.
# This script redoes both, then continues with the originally-planned
# nine-pack pooling-grid completion + one-hop capacity diagnostic.
set -uo pipefail  # no -e: one arm failing must not kill the whole night
cd "$(dirname "$0")/.."
source scripts/activate_data_environment.sh

export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
export PYTHONUNBUFFERED=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

LOG_DIR=scratch/logs
mkdir -p "$LOG_DIR"
REPORT=reports/inspection/stage0_7h_nine_pack_smoke
HARD_STOP_EPOCH=$(date -d '2026-09-08 08:30' +%s)

log() { echo "[$(date -Is)] $*"; }
past_cutoff() { [[ "$(date +%s)" -ge "$HARD_STOP_EPOCH" ]]; }

run_onehop_correctness_smokes() {
  if past_cutoff; then log "past cutoff -- skip one-hop smokes"; return 0; fi
  log "=== redo: one-hop correctness smokes (import bug fixed) ==="
  uv run python -u scripts/run_7h_onehop_correctness_smokes.py --device cuda --epochs 5 \
    || log "WARNING: one-hop correctness smokes failed -- continuing"
}

run_b4_pooling() {
  if past_cutoff; then log "past cutoff -- skip Track B.4"; return 0; fi
  log "=== redo: Track B.4 ATS seed-43 2x2 pooling (never ran -- queue crashed before this) ==="
  bash scripts/run_7h_ats_pooling_s2.sh \
    || log "WARNING: Track B.4 failed -- continuing"
}

train_scalar_pooling_pair() {
  if past_cutoff; then log "past cutoff -- skip scalar pooling pair"; return 0; fi
  log "=== nine-pack parity: scalar mean/max + scalar max/mean (completes 5-combo grid at scale) ==="
  local arms=(
    "configs/experiment/stage0_7h_nine_pack_scalar_mean_max.yaml|stage0-7h-nine-pack-scalar-mean-max|_staging_scalar_mean_max"
    "configs/experiment/stage0_7h_nine_pack_scalar_max_mean.yaml|stage0-7h-nine-pack-scalar-max-mean|_staging_scalar_max_mean"
  )
  local entry cfg run_id staging
  for entry in "${arms[@]}"; do
    if past_cutoff; then log "past cutoff mid-pair -- stop before next arm"; return 0; fi
    IFS='|' read -r cfg run_id staging <<<"$entry"
    log "=== train $run_id ==="
    uv run mbs train cascade \
      --config "$cfg" --run-id "$run_id" --device cuda --skip-if-done \
      --report-dir "$REPORT/$staging" \
      || log "WARNING: $run_id failed -- continuing"
  done
}

train_onehop_wide_diagnostic() {
  if past_cutoff; then log "past cutoff -- skip one-hop wide-capacity diagnostic"; return 0; fi
  log "=== diagnostic: one-hop m-only with rho_hidden_dimension 10->64 (3 folds) ==="
  uv run python - <<'PY' || echo "[queue] wide diagnostic failed -- continuing"
import sys
from pathlib import Path
ROOT = Path(".").resolve()
sys.path.insert(0, str(ROOT / "scripts"))
from run_7g_gene_only_probe import train_flat_region_arm
from mbs.paths import DataPaths

paths = DataPaths.from_environment()
train_flat_region_arm(
    paths=paths,
    config_path=ROOT / "configs/experiment/stage0_7h_nine_pack_m_only_wide.yaml",
    run_prefix="stage0-7h-nine-pack-m-only-wide",
    device="cuda",
    report_dir=ROOT / "reports/inspection/stage0_7h_nine_pack_smoke",
    fold_filter=None,
)
print("[queue] m-only-wide diagnostic done", flush=True)
PY
}

refresh_report() {
  log "=== refresh nine-pack report ==="
  uv run python -u scripts/run_7h_nine_pack_smoke.py --phase report || log "WARNING: report refresh failed"
}

log "=== 7H overnight chain v2 start (CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES, hard stop 08:30 CEST) ==="
run_onehop_correctness_smokes
run_b4_pooling
refresh_report
train_scalar_pooling_pair
refresh_report
train_onehop_wide_diagnostic
refresh_report

log "=== overnight chain done (or hit 08:30 cutoff) ==="
log "GPU 0 idle from here unless manually resumed. Do NOT auto-launch Stage B / Milestone 12 OOF / disease GPU."
log "Reports: $REPORT/analysis.md"
