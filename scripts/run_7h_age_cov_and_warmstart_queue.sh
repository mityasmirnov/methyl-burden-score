#!/usr/bin/env bash
# Milestone 10: (1) real 15-epoch/3-fold age-covariates (tissue+sex-conditioned
# age head) ablation vs the P2-G baseline (13.431 MAE), then (2) two vector
# (region_hidden) warm-start arms testing whether LP-FT from a converged
# scalar_rbs checkpoint closes the gap to P2-G scalar (see the warmstart
# config docstrings for full rationale). Each step is independent; one
# failing must not block the others.
set -uo pipefail
cd "$(dirname "$0")/.."
source scripts/activate_data_environment.sh

export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
export PYTHONUNBUFFERED=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

REPORT=reports/inspection/stage0_7h_nine_pack_smoke
log() { echo "[$(date -Is)] $*"; }

refresh_report() {
  log "=== refresh nine-pack report ==="
  uv run python -u scripts/run_7h_nine_pack_smoke.py --phase report || log "WARNING: report refresh failed"
}

log "=== step 1/3: age-covariates (tissue+sex) ablation, full 15ep/3fold ==="
uv run mbs train cascade \
  --config configs/experiment/stage0_7h_nine_pack_p2_g_age_covariates.yaml \
  --run-id stage0-7h-nine-pack-p2-g-age-covariates --device cuda --skip-if-done \
  --report-dir "$REPORT/_staging_p2_g_age_covariates" \
  || log "WARNING: age-covariates ablation failed -- continuing"
refresh_report

log "=== step 2/3: vector max/max warm-started from P2-G scalar checkpoint ==="
uv run mbs train cascade \
  --config configs/experiment/stage0_7h_nine_pack_vector_max_max_warmstart.yaml \
  --run-id stage0-7h-nine-pack-vector-max-max-warmstart --device cuda --skip-if-done \
  --report-dir "$REPORT/_staging_vector_max_max_warmstart" \
  || log "WARNING: vector max/max warm-start failed -- continuing"
refresh_report

log "=== step 3/3: vector mean/max warm-started from scalar mean/max checkpoint ==="
uv run mbs train cascade \
  --config configs/experiment/stage0_7h_nine_pack_vector_mean_max_warmstart.yaml \
  --run-id stage0-7h-nine-pack-vector-mean-max-warmstart --device cuda --skip-if-done \
  --report-dir "$REPORT/_staging_vector_mean_max_warmstart" \
  || log "WARNING: vector mean/max warm-start failed -- continuing"
refresh_report

log "=== queue done ==="
