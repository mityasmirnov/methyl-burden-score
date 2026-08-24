#!/usr/bin/env bash
# Thin launcher: nohup the 7B convert+progress watcher (idempotent).
# Does not kill an in-flight convert-pack; only replaces the progress watcher.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

LOGDIR="$MBS_ARTIFACT_ROOT/logs/matrix_convert"
mkdir -p "$LOGDIR"
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
NOHUP_LOG="$LOGDIR/7b_hub_matrices_nohup_${STAMP}.log"
PID_FILE="$LOGDIR/7b_hub_matrices_nohup_${STAMP}.pid"

# Stop previous watcher only (never convert-pack / convert_hub_full_packs.sh).
if pgrep -f "scripts/run_7b_hub_matrices_background.sh" >/dev/null 2>&1; then
  echo "Stopping previous 7B watcher (convert-pack left running)…"
  pkill -f "scripts/run_7b_hub_matrices_background.sh" || true
  sleep 1
fi

nohup bash "$REPO_ROOT/scripts/run_7b_hub_matrices_background.sh" >"$NOHUP_LOG" 2>&1 &
PID=$!
printf '%s\n' "$PID" >"$PID_FILE"
ln -sfn "$(basename "$NOHUP_LOG")" "$LOGDIR/7b_hub_matrices_nohup_latest.log"
ln -sfn "$(basename "$PID_FILE")" "$LOGDIR/7b_hub_matrices_nohup_latest.pid"

printf 'Started 7B convert watcher\n'
printf '  pid:     %s\n' "$PID"
printf '  log:     %s\n' "$NOHUP_LOG"
printf '  latest:  %s\n' "$LOGDIR/7b_hub_matrices_nohup_latest.log"
printf '  pidfile: %s\n' "$PID_FILE"
printf 'Track:     bash scripts/status_7b_hub_matrices.sh\n'
printf 'Progress:  reports/inspection/stage0_7b_hub_matrices/progress.md\n'
printf 'Poll:      rg FINALIZE_DONE %s\n' "$LOGDIR/7b_hub_matrices_nohup_latest.log"
