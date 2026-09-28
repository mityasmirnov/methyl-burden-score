#!/usr/bin/env bash
# Remainder after N-light-first in-flight queue: cascade G1 never ran.
# Skip completed N-light train+nested; run cascade → full-width → baseline.
set -uo pipefail
cd "$(dirname "$0")/.."
source scripts/activate_data_environment.sh 2>/dev/null || true

export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export PYTHONUNBUFFERED=1

log() { printf '[%s] %s\n' "$(date -Is)" "$*"; }
run_job() {
  local name="$1"; shift
  log "=== START $name ==="
  if "$@"; then log "=== OK $name ==="; else log "=== FAILED $name (rc=$?) -- continuing ==="; fi
}

log "GPU0 G1 remainder: cascade first (N-light already done), then §3 fill"

run_job "g1-cascade-random-train" \
  uv run python -u scripts/run_12_cascade_oof.py \
    --config configs/experiment/stage0_12b_gene_holdout_cascade_smoke.yaml \
    --run-prefix stage0-12b-gh-cascade-random \
    --device cuda \
    --report-dir reports/inspection/stage0_12b_gene_holdout_cascade \
    --folds 0 --restarts 0 --skip-enet
run_job "g1-cascade-random-nested" \
  uv run python -u scripts/eval_gene_holdout_nested.py \
    --run-id stage0-12b-gh-cascade-random-f0-r0 --fold 0 --force

run_job "full-width-converge" \
  uv run python -u scripts/run_12b_cpgpt_full_width_converge.py
run_job "matched-baseline-6restart" \
  uv run python -u scripts/run_12_nlight_oof.py \
    --config configs/experiment/stage0_7h_nine_pack_m_only_wide.yaml \
    --run-prefix stage0-7h-baseline-multirestart \
    --device cuda \
    --report-dir reports/inspection/stage0_7h_baseline_multirestart \
    --folds 0 --restarts all

log "GPU0 G1 remainder complete"
