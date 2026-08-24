#!/usr/bin/env bash
# Background Milestone 7B convert + progress watcher.
#
# - If convert_hub_full_packs.sh is not already running, start it (skip-if-exists).
# - Poll every INTERVAL_SEC; rewrite progress.md/json + plan Progress block.
# - When all six 7B packs have conversion_stats.json and convert is idle:
#     index-hub-packs, write_stage0_7b_report, catalog refresh-release, final progress.
#
# Usage:
#   bash scripts/convert_hub_full_packs_background.sh   # preferred entry
#   bash scripts/run_7b_hub_matrices_background.sh      # foreground / nohup target
# Track:
#   bash scripts/status_7b_hub_matrices.sh
#   cat reports/inspection/stage0_7b_hub_matrices/progress.md
set -euo pipefail
cd "$(dirname "$0")/.."
# shellcheck disable=SC1091
source scripts/activate_data_environment.sh

export PYTHONUNBUFFERED=1
INTERVAL_SEC="${MBS_7B_PROGRESS_INTERVAL_SEC:-30}"
LOGDIR="${MBS_ARTIFACT_ROOT}/logs/matrix_convert"
mkdir -p "$LOGDIR"
TS="$(date -u +%Y%m%dT%H%M%SZ)"
WATCH_LOG="$LOGDIR/7b_watcher_${TS}.log"
CONVERT_LOG="$LOGDIR/hub_full_packs_${TS}.log"
PID_FILE="$LOGDIR/7b_hub_matrices_${TS}.pid"
LATEST_WATCH="$LOGDIR/7b_watcher_latest.log"
LATEST_PID="$LOGDIR/7b_hub_matrices_latest.pid"
PROGRESS_PY="scripts/update_7b_convert_progress.py"
REPORT_DIR="reports/inspection/stage0_7b_hub_matrices"

echo "$$" >"$PID_FILE"
ln -sfn "$(basename "$PID_FILE")" "$LATEST_PID"
ln -sfn "$(basename "$WATCH_LOG")" "$LATEST_WATCH"
exec > >(tee -a "$WATCH_LOG") 2>&1

echo "=== 7B watcher starting ${TS} pid=$$ interval=${INTERVAL_SEC}s ==="
echo "progress: ${REPORT_DIR}/progress.md"
echo "pidfile: ${PID_FILE}"
echo "latest_log: ${LATEST_WATCH}"

update_progress() {
  uv run python "$PROGRESS_PY" || true
}

convert_running() {
  pgrep -f "scripts/convert_hub_full_packs.sh" >/dev/null 2>&1 \
    || pgrep -f "mbs matrix convert-pack" >/dev/null 2>&1
}

seven_b_done_count() {
  local n=0 family
  for family in ancestry bmi brain blood cancer disease; do
    if [[ -f "data/canonical/matrices/matrix-hub-${family}-full-v1/conversion_stats.json" ]]; then
      n=$((n + 1))
    fi
  done
  echo "$n"
}

# Kick convert if nothing is converting (safe: skip-if-exists).
if ! convert_running; then
  echo "--- no convert process found; launching convert_hub_full_packs.sh"
  nohup bash scripts/convert_hub_full_packs.sh >"$CONVERT_LOG" 2>&1 &
  echo "convert_pid=$! convert_log=$CONVERT_LOG"
  ln -sfn "$(basename "$CONVERT_LOG")" "$LOGDIR/hub_full_packs_latest.log"
else
  echo "--- convert already running; watcher will only track + finalize"
fi

update_progress

finalize() {
  echo "--- finalize: index + inspection report + release refresh"
  uv run mbs matrix index-hub-packs \
    --check-overlap \
    --report-dir "$REPORT_DIR"
  uv run python scripts/write_stage0_7b_report.py
  uv run mbs catalog refresh-release || echo "WARN: catalog refresh-release failed (non-fatal for matrices)"
  update_progress
  echo "FINALIZE_DONE"
}

FINALIZED=0
while true; do
  update_progress
  DONE="$(seven_b_done_count)"
  echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) seven_b_done=${DONE}/6 convert_running=$(convert_running && echo yes || echo no)"
  if [[ "$DONE" -eq 6 ]] && ! convert_running; then
    if [[ "$FINALIZED" -eq 0 ]]; then
      finalize
      FINALIZED=1
    fi
    echo "=== 7B watcher complete (6/6 packs) ==="
    exit 0
  fi
  sleep "$INTERVAL_SEC"
done
