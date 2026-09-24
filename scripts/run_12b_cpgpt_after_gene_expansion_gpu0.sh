#!/usr/bin/env bash
# Takeover waiter: the already-running run_12b_gene_expansion_smoke_gpu0.sh
# (bash PID $WRAPPER_PID) chains straight from the N-light full-gene-panel
# smoke into a P2-G cascade full-panel smoke. Per updated priority
# (2026-09-08 evening), that cascade step is deferred in favor of testing
# CpGPT positional embeddings on N-light first. This script watches the
# N-light training process ($NLIGHT_PID) for exit, kills the wrapper bash
# script before it can launch the cascade step (if the wrapper is still
# alive), then launches the CpGPT smoke on GPU0.
#
# Uses kill -0 on tracked PIDs, not pgrep -f text matching -- a pgrep -f
# self-matching bug in an earlier version of this pattern (see
# milestone-7h-pretrained-mbs-rbs-campaign.md 2026-09-07) wasted >1h of GPU
# idle time by matching the text of later diagnostic commands.
set -uo pipefail
cd "$(dirname "$0")/.."
source scripts/activate_data_environment.sh

NLIGHT_PID="${NLIGHT_PID:?set NLIGHT_PID to the run_12b_gene_expansion_nlight_smoke.py PID}"
WRAPPER_PID="${WRAPPER_PID:?set WRAPPER_PID to the run_12b_gene_expansion_smoke_gpu0.sh bash PID}"
POLL_S="${POLL_S:-30}"

log() { echo "[$(date -Is)] $*"; }

log "watching N-light gene-expansion smoke pid=${NLIGHT_PID} (wrapper=${WRAPPER_PID})"
while kill -0 "$NLIGHT_PID" 2>/dev/null; do
  sleep "$POLL_S"
done
log "N-light gene-expansion smoke pid ${NLIGHT_PID} gone"

if kill -0 "$WRAPPER_PID" 2>/dev/null; then
  log "wrapper ${WRAPPER_PID} still alive -- killing before it launches the cascade full-panel smoke"
  kill "$WRAPPER_PID" 2>/dev/null || true
  sleep 2
  kill -9 "$WRAPPER_PID" 2>/dev/null || true
else
  log "wrapper ${WRAPPER_PID} already gone"
fi

export CUDA_VISIBLE_DEVICES=0
export PYTHONUNBUFFERED=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

log "=== 12b CpGPT smoke: N-light (GPU0) ==="
nvidia-smi --query-gpu=index,memory.used,memory.free --format=csv -i 0 || true
uv run python -u scripts/run_12b_cpgpt_nlight_smoke.py
log "=== CpGPT smoke exit=$? ==="
