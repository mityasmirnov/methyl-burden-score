#!/usr/bin/env bash
# Status for 7G′ optional encoder parity background runs.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

LOGDIR="$MBS_ARTIFACT_ROOT/logs/train"
LATEST_LOG="$LOGDIR/stage0_7g_encoder_parity_latest.log"
LATEST_PID="$LOGDIR/stage0_7g_encoder_parity_latest.pid"
REPORT="${REPORT_DIR:-$REPO_ROOT/reports/inspection/stage0_7g_encoder_parity}"

PREFIXES=(
  stage0-7g-encoder-parity-Flat-G
  stage0-7g-encoder-parity-Hier-G
)

printf '=== 7G′ encoder parity status ===\n'
if [[ -f "$LATEST_PID" ]]; then
  PID="$(cat "$LATEST_PID")"
  printf 'pidfile: %s (pid=%s)\n' "$LATEST_PID" "$PID"
  if kill -0 "$PID" 2>/dev/null; then
    printf 'process: RUNNING\n'
  else
    printf 'process: not running (stale pid)\n'
  fi
else
  printf 'pidfile: (none)\n'
fi
if [[ -f "$LATEST_LOG" ]]; then
  printf 'log: %s\n' "$LATEST_LOG"
  printf '%s\n' '--- last 15 log lines ---'
  tail -n 15 "$LATEST_LOG" || true
else
  printf 'log: (none)\n'
fi
printf '%s\n' '--- fold metrics ---'
for prefix in "${PREFIXES[@]}"; do
  done_n=0
  for i in 0 1 2; do
    if [[ -f "$MBS_ARTIFACT_ROOT/runs/${prefix}-f${i}/metrics.json" ]]; then
      done_n=$((done_n + 1))
    fi
  done
  printf '%s: %s/3 folds\n' "$prefix" "$done_n"
done
if [[ -f "$REPORT/summary.json" ]]; then
  printf 'report: %s\n' "$REPORT"
else
  printf 'report: pending (%s)\n' "$REPORT"
fi
