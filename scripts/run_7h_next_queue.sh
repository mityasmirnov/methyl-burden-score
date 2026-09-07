#!/usr/bin/env bash
# 7H next-job queue (post Track A smoke).
#
# Does NOT interrupt the live nine-pack full refs. Waits for
# scratch/logs/7h_nine_pack_full.pid to exit, then runs Track B.4
# (ATS seed-43 2×2 pooling). Hard-stops before Stage B / OOF / disease GPU.
#
# Usage:
#   nohup bash scripts/run_7h_next_queue.sh >> scratch/logs/7h_next_queue.log 2>&1 &
#   echo $! > scratch/logs/7h_next_queue.pid
set -euo pipefail
cd "$(dirname "$0")/.."
source scripts/activate_data_environment.sh

export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
export PYTHONUNBUFFERED=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

LOG_DIR=scratch/logs
mkdir -p "$LOG_DIR"
FULL_PID_FILE="$LOG_DIR/7h_nine_pack_full.pid"
QUEUE_LOG="$LOG_DIR/7h_next_queue.log"

log() { echo "[$(date -Is)] $*"; }

wait_for_pidfile() {
  local pidf="$1"
  local label="$2"
  if [[ ! -f "$pidf" ]]; then
    log "no $pidf — assuming $label already finished"
    return 0
  fi
  local pid
  pid="$(tr -d '[:space:]' <"$pidf" || true)"
  if [[ -z "$pid" ]]; then
    log "empty $pidf — assuming $label finished"
    return 0
  fi
  if ! kill -0 "$pid" 2>/dev/null; then
    log "$label pid $pid not alive — continuing"
    return 0
  fi
  log "waiting for $label pid=$pid"
  while kill -0 "$pid" 2>/dev/null; do
    sleep 120
  done
  log "$label pid=$pid exited"
}

# Also wait if cascade/full still owns GPU even without pidfile match.
wait_for_ninepack_gpu_owner() {
  while pgrep -af 'run_7h_nine_pack_smoke.py --phase full|stage0-7h-nine-pack-P2-G[^-]|stage0-7h-nine-pack-m-only' \
      | rg -v 'rg |pgrep|next_queue' >/dev/null; do
    log "nine-pack GPU owner still present — sleep 120s"
    sleep 120
  done
}

log "=== 7H next-queue start ==="
wait_for_pidfile "$FULL_PID_FILE" "nine-pack-full"
wait_for_ninepack_gpu_owner

log "=== refresh nine-pack full report ==="
uv run python -u scripts/run_7h_nine_pack_smoke.py --phase report || true

log "=== Track B.4: ATS seed-43 2×2 pooling ==="
bash scripts/run_7h_ats_pooling_s2.sh

log "=== HARD STOP ==="
log "Do not launch Stage B / Milestone 7 OOF / disease-cancer GPU automatically."
log "Next (manual Track C): define disease/cancer case-control labels (false≠control),"
log "  then optional nine-pack trait arm; blood/brain masks need phenotype repair first."
log "Recommended review:"
log "  reports/inspection/stage0_7h_nine_pack_smoke/analysis.md"
log "  reports/inspection/stage0_7h_nine_pack_smoke/trait_adequacy.md"
log "  scratch/logs/7h_ats_pooling_s2.log"
log "=== queue done ==="
