#!/usr/bin/env bash
# Wait for dense-stage1 *S1* to finish on GPU 2, then start N-light 5×6 OOF
# at max VRAM. If the dense queue proceeds to naive max/max transplants
# (step 2/3), stop it — those are not the 10e S1–S4 recipe and GPU 2 is
# now reserved for Milestone 12 N-light.
set -uo pipefail
cd "$(dirname "$0")/.."

LOG="${DENSE_STAGE1_LOG:-scratch/logs/7h_dense_stage1_queue.log}"
QUEUE_PID="${DENSE_STAGE1_QUEUE_PID:-}"
POLL_S="${POLL_S:-30}"

log() { echo "[$(date -Is)] $*"; }

if [[ -z "$QUEUE_PID" ]]; then
  QUEUE_PID=$(pgrep -f 'bash scripts/run_7h_dense_stage1_queue.sh' | head -n 1 || true)
fi

log "watching dense-stage1 queue pid=${QUEUE_PID:-none} log=$LOG"

stop_transplants() {
  log "S1 complete — stopping dense-stage1 queue before max/max transplants"
  if [[ -n "${QUEUE_PID}" ]] && kill -0 "$QUEUE_PID" 2>/dev/null; then
    pkill -P "$QUEUE_PID" 2>/dev/null || true
    kill "$QUEUE_PID" 2>/dev/null || true
    sleep 2
    kill -9 "$QUEUE_PID" 2>/dev/null || true
  fi
  # Any leftover cascade train on GPU 2 from this queue.
  pkill -f 'stage0-7h-nine-pack-scalar-max-max-from-dense-stage1' 2>/dev/null || true
  pkill -f 'stage0-7h-nine-pack-vector-max-max-from-dense-stage1' 2>/dev/null || true
}

while true; do
  if [[ -f "$LOG" ]] && grep -q '=== step 2/3' "$LOG"; then
    stop_transplants
    break
  fi
  if [[ -f "$LOG" ]] && grep -q '=== dense stage-1 queue done ===' "$LOG"; then
    log "dense-stage1 queue finished on its own"
    break
  fi
  if [[ -n "${QUEUE_PID}" ]] && ! kill -0 "$QUEUE_PID" 2>/dev/null; then
    log "dense-stage1 queue pid $QUEUE_PID gone"
    break
  fi
  if [[ -z "${QUEUE_PID}" ]]; then
    if pgrep -f 'stage0-7h-nine-pack-dense-stage1-mean-mean' >/dev/null; then
      :
    else
      log "no dense-stage1 process; GPU 2 assumed free"
      break
    fi
  fi
  sleep "$POLL_S"
done

# Wait until GPU 2 has no remaining dense-stage1 python.
for _ in $(seq 1 60); do
  if pgrep -f 'stage0-7h-nine-pack-dense-stage1-mean-mean' >/dev/null; then
    log "waiting for dense-stage1 python to exit"
    sleep 10
    continue
  fi
  break
done

log "launching N-light OOF on GPU 2"
exec bash scripts/run_12_nlight_oof.sh
