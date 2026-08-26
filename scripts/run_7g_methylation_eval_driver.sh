#!/usr/bin/env bash
# Idempotent 7G driver: cascade folds (skip-if-done) then report writer.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

CONFIG="$REPO_ROOT/configs/experiment/stage0_7g_methylation_eval.yaml"
RUN_ID="stage0-7g-cascade-v1"
REPORT_DIR="$REPO_ROOT/reports/inspection/stage0_7g_methylation_eval"
DEVICE="cuda"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2 ;;
    --run-id) RUN_ID="$2"; shift 2 ;;
    --report-dir) REPORT_DIR="$2"; shift 2 ;;
    --device) DEVICE="$2"; shift 2 ;;
    *) echo "unknown arg: $1" >&2; exit 2 ;;
  esac
done

echo "[7g] cascade train run_id=$RUN_ID device=$DEVICE" >&2
export PYTHONUNBUFFERED=1
uv run mbs train cascade \
  --config "$CONFIG" \
  --run-id "$RUN_ID" \
  --device "$DEVICE" \
  --report-dir "$REPORT_DIR" \
  --skip-if-done

echo "[7g] classical + transparent + report" >&2
uv run python "$REPO_ROOT/scripts/write_7g_methylation_eval_report.py" \
  --config "$CONFIG" \
  --run-id "$RUN_ID" \
  --report-dir "$REPORT_DIR"

echo "[7g] done report=$REPORT_DIR" >&2
