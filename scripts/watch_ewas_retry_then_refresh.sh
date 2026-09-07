#!/usr/bin/env bash
# When clean EWAS retry finishes: resummarize failures (+ catalog refresh if needed).
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1091
source "$ROOT/scripts/activate_data_environment.sh"
cd "$ROOT"
OUTER="$ROOT/scratch/logs/ewas_db_retry_clean.outer.log"
WATCHLOG="$ROOT/scratch/logs/ewas_retry_watch.log"
mkdir -p "$ROOT/scratch/logs"
log() { printf '%s %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a "$WATCHLOG"; }

log "watching for retry_ewas_db_download_failures.sh exit"
while pgrep -f 'scripts/retry_ewas_db_download_failures.sh' >/dev/null 2>&1; do
  sleep 60
done
log "retry process gone"
sleep 5
log "summarize failures"
uv run python scripts/summarize_ewas_db_download_failures.py | tee -a "$WATCHLOG"
if [[ -f "$OUTER" ]] && grep -q 'catalog-refresh-release\|post_ewas_datahub' "$OUTER" 2>/dev/null; then
  log "post-hook likely ran with retry; skip duplicate refresh"
else
  log "running catalog-refresh-release (MBS_SKIP_ATLAS_SEED=1)"
  MBS_SKIP_ATLAS_SEED=1 make catalog-refresh-release | tee -a "$WATCHLOG"
fi
log "watch done"
