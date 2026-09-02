#!/usr/bin/env bash
# Status for 7G′ Stage A background gene-only probe (pid, log tail, arm progress).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

LOGDIR="$MBS_ARTIFACT_ROOT/logs/train"
LATEST_LOG="$LOGDIR/stage0_7g_gene_only_latest.log"
LATEST_PID="$LOGDIR/stage0_7g_gene_only_latest.pid"
REPORT="${REPORT_DIR:-$REPO_ROOT/reports/inspection/stage0_7g_gene_only_probe}"

ARMS=(
  stage0-7g-gene-probe-P2-G
  stage0-7g-gene-probe-P4-G
  stage0-7g-gene-probe-P5-G-max
  stage0-7g-gene-probe-P5-G-mean
  stage0-7g-gene-probe-P2-orphan-ablation
)

printf '=== 7G′ Stage A gene-only probe status ===\n'
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
RUNNING="$(pgrep -af 'run_7g_gene_only_probe.py' 2>/dev/null | grep -v status_7g || true)"
if [[ -n "$RUNNING" ]]; then
  printf '%s\n' '--- active arm jobs ---'
  printf '%s\n' "$RUNNING"
fi
CLASSICAL_LOG="$REPO_ROOT/scratch/logs/gene_probe_classical_g.log"
if [[ -f "$CLASSICAL_LOG" ]]; then
  printf 'classical log: %s\n' "$CLASSICAL_LOG"
  tail -n 3 "$CLASSICAL_LOG" 2>/dev/null || true
fi
if [[ -f "$REPORT/per_arm/C-mvalue-classical-G.json" ]]; then
  printf 'classical per_arm: present\n'
else
  printf 'classical per_arm: pending\n'
fi
if [[ -f "$LATEST_LOG" ]]; then
  printf 'log: %s\n' "$LATEST_LOG"
  printf '%s\n' '--- last 20 log lines ---'
  tail -n 20 "$LATEST_LOG" || true
else
  printf 'log: (none)\n'
fi
printf '%s\n' '--- cascade run fold metrics ---'
for run_id in "${ARMS[@]}"; do
  RUN_ROOT="$MBS_ARTIFACT_ROOT/runs/$run_id"
  done_n=0
  for i in 0 1 2; do
    if [[ -f "$RUN_ROOT/fold_$i/metrics.json" ]]; then
      done_n=$((done_n + 1))
    fi
  done
  printf '%s: %s/3 folds\n' "$run_id" "$done_n"
done
if [[ -f "$REPORT/analysis.md" ]]; then
  printf 'report: %s (analysis.md present)\n' "$REPORT"
elif [[ -f "$REPORT/summary.json" ]]; then
  printf 'report: %s (summary only; analysis pending)\n' "$REPORT"
else
  printf 'report: pending (%s)\n' "$REPORT"
fi
