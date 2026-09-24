#!/usr/bin/env bash
# Wait for Milestone 12 N-light OOF (GPU 2) to finish, then launch the
# P2-G cascade 5x6 OOF on GPU 2. Uses kill -0 on the tracked PID, not
# pgrep -f text matching -- a pgrep -f self-matching bug in an earlier
# version of this pattern (see milestone-7h-pretrained-mbs-rbs-campaign.md
# 2026-09-07) wasted >1h of GPU idle time by matching the text of later
# diagnostic commands that merely mentioned the watched script's filename.
set -uo pipefail
cd "$(dirname "$0")/.."
source scripts/activate_data_environment.sh

NLIGHT_PID="${NLIGHT_PID:?set NLIGHT_PID to the run_12_nlight_oof.sh bash PID}"
POLL_S="${POLL_S:-60}"

log() { echo "[$(date -Is)] $*"; }

log "watching N-light OOF pid=${NLIGHT_PID}"
while kill -0 "$NLIGHT_PID" 2>/dev/null; do
  sleep "$POLL_S"
done
log "N-light OOF pid ${NLIGHT_PID} gone -- launching cascade OOF on GPU 2"

export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-2}"
export PYTHONUNBUFFERED=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

REPORT=reports/inspection/stage0_12_cascade_oof
mkdir -p "$REPORT"
log "=== Milestone 12 cascade OOF GPU=${CUDA_VISIBLE_DEVICES} ==="
nvidia-smi --query-gpu=index,name,memory.total,memory.used,utilization.gpu --format=csv || true

uv run python -u scripts/run_12_cascade_oof.py \
  --config configs/experiment/stage0_12_cascade_oof.yaml \
  --run-prefix stage0-12-cascade-oof \
  --device cuda \
  --report-dir "$REPORT"

log "=== cascade OOF runner exited $? ==="
