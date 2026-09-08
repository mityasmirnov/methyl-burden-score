#!/usr/bin/env bash
# Milestone 10: dense-gradient stage-1 RBS pretrain + two stage-2 transplants
# (scalar max/max, vector max/max), on GPU 2. See config docstrings for full
# rationale: region_rho never gets a dense gradient under max-pooling, and
# RBS-based classical probes already beat mbs_e2e everywhere -- this tests
# whether a more thoroughly-trained encoder pushes either further.
set -uo pipefail
cd "$(dirname "$0")/.."
source scripts/activate_data_environment.sh

export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-2}"
export PYTHONUNBUFFERED=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

REPORT=reports/inspection/stage0_7h_nine_pack_smoke
log() { echo "[$(date -Is)] $*"; }

refresh_report() {
  log "=== refresh nine-pack report ==="
  uv run python -u scripts/run_7h_nine_pack_smoke.py --phase report || log "WARNING: report refresh failed"
}

log "=== step 1/3: dense-gradient stage-1 pretrain (cpg_pool: mean, region_pool: mean) ==="
uv run mbs train cascade \
  --config configs/experiment/stage0_7h_nine_pack_dense_stage1_mean_mean.yaml \
  --run-id stage0-7h-nine-pack-dense-stage1-mean-mean --device cuda --skip-if-done \
  --report-dir "$REPORT/_staging_dense_stage1_mean_mean" \
  || { log "FATAL: dense stage-1 pretrain failed -- stage-2 transplants need this checkpoint, stopping queue"; exit 1; }

log "=== step 2/3: scalar max/max warm-started from dense stage-1 ==="
uv run mbs train cascade \
  --config configs/experiment/stage0_7h_nine_pack_scalar_max_max_from_dense_stage1.yaml \
  --run-id stage0-7h-nine-pack-scalar-max-max-from-dense-stage1 --device cuda --skip-if-done \
  --report-dir "$REPORT/_staging_scalar_max_max_from_dense_stage1" \
  || log "WARNING: scalar from-dense-stage1 failed -- continuing"
refresh_report

log "=== step 3/3: vector max/max warm-started from dense stage-1 ==="
uv run mbs train cascade \
  --config configs/experiment/stage0_7h_nine_pack_vector_max_max_from_dense_stage1.yaml \
  --run-id stage0-7h-nine-pack-vector-max-max-from-dense-stage1 --device cuda --skip-if-done \
  --report-dir "$REPORT/_staging_vector_max_max_from_dense_stage1" \
  || log "WARNING: vector from-dense-stage1 failed -- continuing"
refresh_report

log "=== dense stage-1 queue done ==="
