#!/usr/bin/env bash
# Status for 7G′ Stage B matched-panel benchmark.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

LOGDIR="$MBS_ARTIFACT_ROOT/logs/train"
LATEST_LOG="$LOGDIR/stage0_7g_prime_stage_b_latest.log"
LATEST_PID="$LOGDIR/stage0_7g_prime_stage_b_latest.pid"
REPORT="${REPORT_DIR:-$REPO_ROOT/reports/inspection/stage0_7g_prime_matched_probe}"

printf '=== 7G′ Stage B matched-panel status ===\n'
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
if [[ -f "$REPORT/per_arm/C-mvalue-enetS.json" ]]; then
  printf 'C-mvalue-enetS: done\n'
else
  printf 'C-mvalue-enetS: pending\n'
fi
for i in 0 1 2; do
  cascade_done=
  light_done=
  if [[ -d "$REPORT/_staging_N_cascade_S_fold_${i}" ]]; then
    cascade_done=yes
  fi
  if [[ -d "$REPORT/_staging_N_light_type_fold_${i}" ]]; then
    light_done=yes
  fi
  printf 'fold_%s: N-cascade-S=%s N-light-type=%s\n' "$i" "${cascade_done:-pending}" "${light_done:-pending}"
done
if [[ -f "$REPORT/analysis.md" ]]; then
  printf 'report: %s (analysis.md present)\n' "$REPORT"
elif [[ -f "$REPORT/summary.json" ]]; then
  printf 'report: %s (summary only)\n' "$REPORT"
else
  printf 'report: pending (%s)\n' "$REPORT"
fi
