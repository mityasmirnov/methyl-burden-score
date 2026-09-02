#!/usr/bin/env bash
# Post-download hook: failure summary + deepmat-data-v1 catalog refresh.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

LOG="${EWAS_DATAHUB_LOG:-$MBS_ARTIFACT_ROOT/logs/downloads/ewas_datahub_EWAS_db.log}"

printf '=== post EWAS DataHub download hook ===\n'
printf 'log=%s\n' "$LOG"

if [[ -f "$LOG" ]]; then
  uv run python "$SCRIPT_DIR/summarize_ewas_db_download_failures.py" --log "$LOG" || {
    printf 'WARN: failure summary failed (continuing to catalog refresh)\n' >&2
  }
else
  printf 'WARN: download log missing (%s); skipping failure summary\n' "$LOG" >&2
fi

printf '=== catalog refresh-release ===\n'
make catalog-refresh-release

printf '=== post hook complete ===\n'
