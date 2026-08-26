#!/usr/bin/env bash
# Status for Milestone 7G background eval (pid, log tail, fold skip-if-done).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

RUN_ID="${RUN_ID:-stage0-7g-cascade-v1}"
LOGDIR="$MBS_ARTIFACT_ROOT/logs/train"
LATEST_LOG="$LOGDIR/stage0_7g_latest.log"
LATEST_PID="$LOGDIR/stage0_7g_latest.pid"
RUN_ROOT="$MBS_ARTIFACT_ROOT/runs/$RUN_ID"
REPORT="$REPO_ROOT/reports/inspection/stage0_7g_methylation_eval"

printf '=== Milestone 7G status ===\n'
printf 'run_id: %s\n' "$RUN_ID"
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
  printf '%s\n' '--- last 20 log lines ---'
  tail -n 20 "$LATEST_LOG" || true
else
  printf 'log: (none)\n'
fi
printf '%s\n' '--- fold score manifests ---'
for i in 0 1 2; do
  MAN="$RUN_ROOT/fold_$i/scores/score_manifest.json"
  MET="$RUN_ROOT/fold_$i/metrics.json"
  if [[ -f "$MAN" ]]; then
    printf 'fold_%s: DONE (%s)\n' "$i" "$MAN"
  elif [[ -f "$MET" ]]; then
    printf 'fold_%s: metrics only (%s)\n' "$i" "$MET"
  else
    printf 'fold_%s: pending\n' "$i"
  fi
done
if [[ -f "$REPORT/analysis.md" ]]; then
  printf 'report: %s (analysis.md present)\n' "$REPORT"
else
  printf 'report: pending (%s)\n' "$REPORT"
fi
