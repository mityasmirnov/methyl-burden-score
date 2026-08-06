#!/usr/bin/env bash
# Launch full Hub pack converts in the background via nohup (Milestone 5d).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

LOGDIR="$MBS_ARTIFACT_ROOT/logs/matrix_convert"
mkdir -p "$LOGDIR"
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
NOHUP_LOG="$LOGDIR/hub_full_packs_nohup_${STAMP}.log"
PID_FILE="$LOGDIR/hub_full_packs_nohup_${STAMP}.pid"

nohup bash "$REPO_ROOT/scripts/convert_hub_full_packs.sh" >"$NOHUP_LOG" 2>&1 &
PID=$!
printf '%s\n' "$PID" >"$PID_FILE"

printf 'Started background Hub full-pack convert\n'
printf '  pid:  %s\n' "$PID"
printf '  log:  %s\n' "$NOHUP_LOG"
printf '  pidfile: %s\n' "$PID_FILE"
printf 'Poll: rg CONVERT_DONE %s\n' "$NOHUP_LOG"
