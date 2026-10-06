#!/usr/bin/env bash
# GPU0 queue — TODO §3.5 within-gene sampler smoke (unblocks full-width).
#
# Do NOT re-queue dense full-width converge (killed 23h52m / ep8 / ~90h projected).
# GPU0 single-owner: pgrep -af run_gpu0_queue before starting.
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

if pgrep -af 'run_gpu0_queue' | grep -v "$$" | grep -v grep >/dev/null; then
  log "ABORT: another run_gpu0_queue is already running"
  pgrep -af 'run_gpu0_queue' | grep -v grep || true
  exit 1
fi

log "GPU0 §3.5 sampler smoke on device ${CUDA_VISIBLE_DEVICES}"

run_job "s35-sampler-smoke" \
  uv run python -u scripts/run_12b_sampler_smoke.py

log "GPU0 §3.5 queue complete"
