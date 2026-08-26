#!/usr/bin/env bash
# Background Milestone 7G eval: cascade (3 folds, skip-if-done) → classical/transparent → report.
# Survives Cursor/session disconnect via nohup. Re-run to resume completed folds.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

RUN_ID="${RUN_ID:-stage0-7g-cascade-v1}"
CONFIG="${CONFIG:-$REPO_ROOT/configs/experiment/stage0_7g_methylation_eval.yaml}"
REPORT_DIR="${REPORT_DIR:-$REPO_ROOT/reports/inspection/stage0_7g_methylation_eval}"
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
DEVICE="${DEVICE:-cuda}"

LOGDIR="$MBS_ARTIFACT_ROOT/logs/train"
mkdir -p "$LOGDIR"
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
NOHUP_LOG="$LOGDIR/stage0_7g_${STAMP}.log"
PID_FILE="$LOGDIR/stage0_7g_${STAMP}.pid"
LATEST_LOG="$LOGDIR/stage0_7g_latest.log"
LATEST_PID="$LOGDIR/stage0_7g_latest.pid"

DRIVER="$REPO_ROOT/scripts/run_7g_methylation_eval_driver.sh"

nohup bash "$DRIVER" \
  --config "$CONFIG" \
  --run-id "$RUN_ID" \
  --report-dir "$REPORT_DIR" \
  --device "$DEVICE" \
  >"$NOHUP_LOG" 2>&1 &
PID=$!
printf '%s\n' "$PID" >"$PID_FILE"
ln -sfn "$NOHUP_LOG" "$LATEST_LOG"
ln -sfn "$PID_FILE" "$LATEST_PID"

printf 'Started background Milestone 7G methylation eval\n'
printf '  run_id: %s\n' "$RUN_ID"
printf '  pid:    %s\n' "$PID"
printf '  log:    %s\n' "$NOHUP_LOG"
printf '  pidfile:%s\n' "$PID_FILE"
printf '  GPU:    CUDA_VISIBLE_DEVICES=%s\n' "$CUDA_VISIBLE_DEVICES"
printf 'Status: bash scripts/status_7g_methylation_eval.sh\n'
printf 'Poll:   tail -f %s\n' "$LATEST_LOG"
