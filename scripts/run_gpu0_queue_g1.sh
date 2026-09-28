#!/usr/bin/env bash
# GPU0 queue — GATE G1 gene-holdout first, per TODO_PIPELINE section 2
# ("NEXT — GATE G1 gene-holdout (do this now)").
#
# G1 is a GATE: if `φ`/`ρ` do not transfer to genes they never trained on, the
# whole gene-agnostic premise (and the full-width/product path built on it) is
# in question. So gating work runs before optimisation work — the previous queue's
# full-width converge and matched baseline are TODO 3.3 / 3.4 and are re-queued
# at the end here.
#
# Order note: TODO lists cascade (2.1) before N-light (2.2); this runs **N-light
# first**, per the user's explicit choice when asked ("N-light first, cascade
# after") — cheaper, and it establishes gene-agnosticism before spending on the
# slower product arm. Flagged rather than silently reordered.
#
# Each arm: train on ~80% of genes, score the disjoint ~20%, then post-hoc
# nested enet on heldout MBS columns the encoder never saw.
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

log "GPU0 G1 queue starting on device ${CUDA_VISIBLE_DEVICES}"

# --- 2.2 N-light gene-holdout (random split) ---
# NOTE: run_12b_gene_expansion_nlight_smoke.py has NO argparse -- passing
# --config to it is silently ignored and trains the wrong config. Use a
# dedicated driver with the config baked in.
run_job "g1-nlight-random-train" \
  uv run python -u scripts/run_12b_gh_nlight_random.py
run_job "g1-nlight-random-nested" \
  uv run python -u scripts/eval_gene_holdout_nested.py \
    --run-id stage0-12b-gh-nlight-random-f0 --force

# --- 2.1 Cascade gene-holdout (random split) ---
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

# --- 3.3 / 3.4 deprioritised work, so the card never idles ---
run_job "full-width-converge" \
  uv run python -u scripts/run_12b_cpgpt_full_width_converge.py
run_job "matched-baseline-6restart" \
  uv run python -u scripts/run_12_nlight_oof.py \
    --config configs/experiment/stage0_7h_nine_pack_m_only_wide.yaml \
    --run-prefix stage0-7h-baseline-multirestart \
    --device cuda \
    --report-dir reports/inspection/stage0_7h_baseline_multirestart \
    --folds 0 --restarts all

log "GPU0 G1 queue complete"
