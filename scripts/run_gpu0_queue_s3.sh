#!/usr/bin/env bash
# GPU0 queue — TODO_PIPELINE section 3 ("AFTER G1 PASSES — lock product path").
#
# G1 random holdout is DONE on both arms and PASSES once the feature-count
# confound is removed (see the matched-n control in gene_holdout_eval.json:
# N-light matched train 11.539 vs heldout 12.047 age MAE, i.e. +0.51 not +3.27).
# So section 2 is closed and this queue runs section 3.
#
# NOT queued, deliberately:
#   3.3 full-width converge — KILLED at 23h52m / epoch 8 of 30. It was still
#       overfitting (val_loss 13.9 -> 42.2) despite lr 5e-4 / patience 10 /
#       dropout 0.3 / wd 1e-3, and projected ~90h at ~3h/epoch (batch 55 forced
#       by VRAM = 622 steps/epoch). Re-queue only after 3.5 (within-gene CpG
#       sampler) cuts edges per sample so the batch can grow. Running it again
#       as-is just burns the card.
#
# GPU0 single-owner rule: only one queue script at a time. On 2026-09-28 two
# queues overlapped on this card (mine plus run_gpu0_queue_g1_remainder.sh).
# Check `pgrep -af run_gpu0_queue` before starting another.
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

log "GPU0 section-3 queue starting on device ${CUDA_VISIBLE_DEVICES}"

# --- 3.2 DeepRVAT-aligned seed gene set (product train-set direction) ---
run_job "s3-nlight-seed-train" \
  uv run python -u scripts/run_12b_gh_nlight_seed.py
run_job "s3-nlight-seed-nested" \
  uv run python -u scripts/eval_gene_holdout_nested.py \
    --run-id stage0-12b-gh-nlight-seed-f0 --force

run_job "s3-cascade-seed-train" \
  uv run python -u scripts/run_12_cascade_oof.py \
    --config configs/experiment/stage0_12b_gene_holdout_cascade_seed_smoke.yaml \
    --run-prefix stage0-12b-gh-cascade-seed \
    --device cuda \
    --report-dir reports/inspection/stage0_12b_gene_holdout_cascade \
    --folds 0 --restarts 0 --skip-enet
run_job "s3-cascade-seed-nested" \
  uv run python -u scripts/eval_gene_holdout_nested.py \
    --run-id stage0-12b-gh-cascade-seed-f0-r0 --fold 0 --force

# --- 3.4 Matched CpGPT-off baseline x6 (rigour: CpGPT effect size is vs n=1) ---
run_job "s3-matched-baseline-6restart" \
  uv run python -u scripts/run_12_nlight_oof.py \
    --config configs/experiment/stage0_7h_nine_pack_m_only_wide.yaml \
    --run-prefix stage0-7h-baseline-multirestart \
    --device cuda \
    --report-dir reports/inspection/stage0_7h_baseline_multirestart \
    --folds 0 --restarts all

log "GPU0 section-3 queue complete"
