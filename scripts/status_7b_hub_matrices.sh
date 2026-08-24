#!/usr/bin/env bash
# One-shot Milestone 7B convert status (updates progress docs, prints summary).
set -euo pipefail
cd "$(dirname "$0")/.."
# shellcheck disable=SC1091
source scripts/activate_data_environment.sh
uv run python scripts/update_7b_convert_progress.py >/dev/null
PROGRESS="reports/inspection/stage0_7b_hub_matrices/progress.md"
echo "=== 7B Hub matrices progress ==="
cat "$PROGRESS"
echo
echo "Convert processes:"
pgrep -af "mbs matrix convert-pack|convert_hub_full_packs" | grep -v pgrep || echo "  (none)"
echo
echo "Watcher:"
pgrep -af "run_7b_hub_matrices_background" | grep -v pgrep || echo "  (none — start with: bash scripts/convert_hub_full_packs_background.sh)"
echo
LATEST_LOG="${MBS_ARTIFACT_ROOT}/logs/matrix_convert/7b_watcher_latest.log"
if [[ -L "$LATEST_LOG" || -f "$LATEST_LOG" ]]; then
  echo "Watcher log (tail): $LATEST_LOG"
  tail -n 8 "$LATEST_LOG" || true
fi
