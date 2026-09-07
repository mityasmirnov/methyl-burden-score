#!/usr/bin/env bash
# Milestone 10 overnight GPU-0 extension (2026-09-07 -> 2026-09-08 09:30 CEST).
#
# User directive: keep GPU 0 running at max capacity nonstop until 09:30
# tomorrow morning; if a decision is needed, make it and report afterward
# rather than blocking on confirmation.
#
# Waits for the already-running scripts/run_7h_next_queue.sh chain (vector
# RBS -> one-hop correctness smokes -> Track B.4) to finish, then:
#   1) Completes nine-pack parity with the ATS-scale 5-combo pooling grid
#      (P2-G scalar max/max + vector mean/max + vector max/max already done
#      or in flight; scalar mean/max + scalar max/mean are the two missing
#      combos -- this is the honest "does pooling choice matter at scale"
#      answer Milestone 11 is gated on, not a new/speculative direction).
#   2) Diagnoses the m-only-vs-P2-G gap that widened at nine-pack scale
#      (0.273 vs 0.355 F1, vs near-parity 0.372/0.379 on ATS): reruns
#      one-hop m-only with rho_hidden_dimension 10->64 (the pooled-
#      representation bottleneck was the narrowest layer in the network by
#      a wide margin -- a concrete, single-variable capacity hypothesis).
# Hard stop at 08:30 CEST (60 min buffer before 09:30) -- does not launch
# any new job past that point, but does not kill anything already running.
# Never launches Stage B / Milestone 12 OOF / disease-cancer GPU -- those
# guardrails hold regardless of time remaining.
set -euo pipefail
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

wait_for_next_queue() {
  log "waiting for scripts/run_7h_next_queue.sh to finish (vector RBS -> one-hop smokes -> B.4)"
  while pgrep -f 'run_7h_next_queue.sh|run_7h_queue_handoff.sh' | grep -v $$ >/dev/null 2>&1; do
    sleep 60
  done
  # Also wait for whatever GPU training process it left running.
  while pgrep -af 'mbs train cascade|run_7h_onehop_correctness_smokes|run_7h_ats_pooling_s2' \
      | grep -v -e "$0" -e 'pgrep' >/dev/null 2>&1; do
    sleep 60
  done
  log "run_7h_next_queue.sh chain finished -- GPU 0 should be free"
}

train_scalar_pooling_pair() {
  if past_cutoff; then
    log "past 08:30 cutoff -- skipping scalar pooling pair"
    return 0
  fi
  log "=== nine-pack parity: scalar mean/max + scalar max/mean (completes 5-combo grid at scale) ==="
  local arms=(
    "configs/experiment/stage0_7h_nine_pack_scalar_mean_max.yaml|stage0-7h-nine-pack-scalar-mean-max|_staging_scalar_mean_max"
    "configs/experiment/stage0_7h_nine_pack_scalar_max_mean.yaml|stage0-7h-nine-pack-scalar-max-mean|_staging_scalar_max_mean"
  )
  local entry cfg run_id staging
  for entry in "${arms[@]}"; do
    if past_cutoff; then
      log "past 08:30 cutoff mid-pair -- stopping before next arm"
      return 0
    fi
    IFS='|' read -r cfg run_id staging <<<"$entry"
    log "=== train $run_id ==="
    uv run mbs train cascade \
      --config "$cfg" \
      --run-id "$run_id" \
      --device cuda \
      --skip-if-done \
      --report-dir "$REPORT/$staging" \
      || log "WARNING: $run_id failed -- continuing"
  done
}

train_onehop_wide_diagnostic() {
  if past_cutoff; then
    log "past 08:30 cutoff -- skipping one-hop wide-capacity diagnostic"
    return 0
  fi
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

log "=== 7H overnight extension start (CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES, hard stop 08:30 CEST) ==="
wait_for_next_queue
refresh_report
train_scalar_pooling_pair
refresh_report
train_onehop_wide_diagnostic
refresh_report

log "=== overnight extension done (or hit 08:30 cutoff) ==="
log "GPU 0 idle from here unless manually resumed. Do NOT auto-launch Stage B / Milestone 12 OOF / disease GPU."
log "Reports: $REPORT/analysis.md"
