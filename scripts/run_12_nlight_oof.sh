#!/usr/bin/env bash
# Milestone 12 N-light 5×6 OOF on GPU 2 (max VRAM). Exclusive: do not share
# this GPU. Recipe: 30 epochs, patience 15, n>200 disease/cancer aux heads.
set -uo pipefail
cd "$(dirname "$0")/.."
source scripts/activate_data_environment.sh

export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-2}"
export PYTHONUNBUFFERED=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
# Exclusive GPU2 — do not leave headroom for a sibling trainer.
export MBS_CASCADE_GPU_SHARE=1
export MBS_CASCADE_GPU_RESERVED_MIB=512

REPORT=reports/inspection/stage0_12_nlight_oof
mkdir -p "$REPORT"
log() { echo "[$(date -Is)] $*"; }

log "=== Milestone 12 N-light OOF GPU=${CUDA_VISIBLE_DEVICES} ==="
nvidia-smi --query-gpu=index,name,memory.total,memory.used,utilization.gpu --format=csv || true

uv run python -u scripts/run_12_nlight_oof.py \
  --config configs/experiment/stage0_12_nlight_oof.yaml \
  --run-prefix stage0-12-nlight-oof \
  --device cuda \
  --report-dir "$REPORT" \
  "$@"

log "=== N-light OOF runner exited $? ==="
