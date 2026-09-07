#!/usr/bin/env bash
# Handoff: wait for in-flight m-only-f0 repair to finish, then replace the
# old keeper process with the updated priority queue (without killing the
# active flat train).
set -euo pipefail
cd "$(dirname "$0")/.."
LOG=scratch/logs/7h_queue_handoff.log
mkdir -p scratch/logs
exec >>"$LOG" 2>&1
log() { echo "[$(date -Is)] $*"; }

METRICS=artifacts/runs/stage0-7h-nine-pack-m-only-f0/metrics.json
OLD_PID_FILE=scratch/logs/7h_next_queue.pid

log "handoff waiting for full-budget m-only-f0 metrics"
while true; do
  if [[ -f "$METRICS" ]]; then
    if python3 - <<'PY'
import json
from pathlib import Path
b=json.loads(Path("artifacts/runs/stage0-7h-nine-pack-m-only-f0/metrics.json").read_text())
raise SystemExit(0 if len(b.get("history") or []) > 3 else 1)
PY
    then
      log "m-only-f0 full-budget metrics present"
      break
    fi
  fi
  sleep 30
done

# If old queue still alive and not yet past vector, stop it so new script runs.
if [[ -f "$OLD_PID_FILE" ]]; then
  old=$(tr -d '[:space:]' <"$OLD_PID_FILE" || true)
  if [[ -n "$old" ]] && kill -0 "$old" 2>/dev/null; then
    cmd=$(ps -p "$old" -o cmd= || true)
    if echo "$cmd" | rg -q 'run_7h_next_queue'; then
      log "stopping old keeper pid=$old after m-only-f0 done"
      kill "$old" || true
      # wait for child uv/python of that queue to exit if any leftover report step
      sleep 5
    fi
  fi
fi

# Do not start if vector already running under another owner
if pgrep -af 'stage0-7h-nine-pack-vector|run_7h_onehop_correctness|run_7h_ats_pooling_s2' \
    | rg -v 'rg |pgrep|handoff' >/dev/null; then
  log "downstream GPU work already running — not relaunching"
  exit 0
fi

nohup bash scripts/run_7h_next_queue.sh >> scratch/logs/7h_next_queue.log 2>&1 &
echo $! > "$OLD_PID_FILE"
log "relaunched keeper pid=$(cat "$OLD_PID_FILE")"
