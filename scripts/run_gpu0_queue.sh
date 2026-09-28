#!/usr/bin/env bash
# GPU0 work queue: run jobs back-to-back so the card never idles.
#
# Ordering rationale (not TODO numbering) -- the goal is the best model, so
# model-advancing jobs come before rigour-only jobs:
#   1. Cascade CpGPT smoke   -- short, and it unblocks cascade OOF (the product
#                               campaign). Plumbing landed but is unproven at
#                               scale, so this is the cheapest way to de-risk it.
#   2. Full-width converge   -- the 19,554-gene run had the best tissue (0.388)
#                               and sex (0.962) on record but never converged
#                               (early-stopped ep 9/16, best ep 4). Most likely
#                               single source of a better model.
#   3. Matched baseline x6   -- CpGPT off, 65k. Rigour: the CpGPT effect size
#                               currently rests on an n=1 baseline. Does not
#                               improve the model, so it runs last.
#
# Each job is independent; a failure is logged and the queue continues to the
# next one rather than losing the remaining GPU time.
set -uo pipefail
cd "$(dirname "$0")/.."
source scripts/activate_data_environment.sh 2>/dev/null || true

export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True   # project convention
export PYTHONUNBUFFERED=1

log() { printf '[%s] %s\n' "$(date -Is)" "$*"; }

run_job() {
  local name="$1"; shift
  log "=== START $name ==="
  if "$@"; then
    log "=== OK $name ==="
  else
    log "=== FAILED $name (rc=$?) -- continuing to next job ==="
  fi
}

log "GPU0 queue starting on device ${CUDA_VISIBLE_DEVICES}"
nvidia-smi --query-gpu=index,name,memory.used,memory.total --format=csv || true

# 1. Cascade + CpGPT plumbing smoke (fold 0, 6 epochs) -- proves the new
#    cascade static-feature path trains and exports scores at real scale.
run_job "cascade-cpgpt-smoke" \
  uv run python -u scripts/run_12_cascade_oof.py \
    --config configs/experiment/stage0_12_cascade_cpgpt_smoke.yaml \
    --run-prefix stage0-12-cascade-cpgpt-smoke \
    --device cuda \
    --report-dir reports/inspection/stage0_12_cascade_cpgpt_smoke \
    --folds 0 --restarts 0

# 2. Full-width (374k col / 19,554 gene) convergence attempt.
run_job "full-width-converge" \
  uv run python -u scripts/run_12b_cpgpt_full_width_converge.py

# 3. Matched no-CpGPT baseline, 6 restarts on fold 0 -- same runner/protocol as
#    the CpGPT multirestart, so the two distributions are directly comparable.
run_job "matched-baseline-6restart" \
  uv run python -u scripts/run_12_nlight_oof.py \
    --config configs/experiment/stage0_7h_nine_pack_m_only_wide.yaml \
    --run-prefix stage0-7h-baseline-multirestart \
    --device cuda \
    --report-dir reports/inspection/stage0_7h_baseline_multirestart \
    --folds 0 --restarts all

log "GPU0 queue complete"
