#!/usr/bin/env bash
# Retry EWAS_db GSM files listed in ewas_db_retry_manifest.tsv (missing after wget WARN).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

MANIFEST="${1:-$MBS_ARTIFACT_ROOT/logs/downloads/ewas_db_retry_manifest.tsv}"
DEST="$MBS_DATA_ROOT/raw/ewas_datahub/EWAS_db"
HTTP_ROOT="https://download.cncb.ac.cn/ewas/datahub/EWAS_db"
LOGDIR="$MBS_ARTIFACT_ROOT/logs/downloads"
mkdir -p "$LOGDIR"
LOG="${EWAS_DATAHUB_RETRY_LOG:-$LOGDIR/ewas_db_retry_$(date -u +%Y%m%dT%H%M%SZ).log}"

if [[ ! -f "$MANIFEST" ]]; then
  printf 'Manifest not found: %s\nRun: uv run python scripts/summarize_ewas_db_download_failures.py\n' "$MANIFEST" >&2
  exit 1
fi

log() { printf '%s %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a "$LOG"; }

log "retry manifest=$MANIFEST dest=$DEST"
tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r study file; do
  [[ -z "$study" || -z "$file" ]] && continue
  mkdir -p "$DEST/$study"
  out="$DEST/$study/$file"
  if [[ -f "$out" && -s "$out" ]]; then
    continue
  fi
  url="${HTTP_ROOT}/${study}/${file}"
  log "wget $study/$file"
  wget -c --tries=3 --retry-connrefused --waitretry=10 --timeout=60 \
    -q -O "$out" "$url" >>"$LOG" 2>&1 || {
    printf 'WARN: failed %s/%s\n' "$study" "$file" | tee -a "$LOG" >&2
  }
done

log "retry pass complete; run post hook or summarize script to refresh manifest"
if [[ "${EWAS_DATAHUB_SKIP_POST_HOOK:-0}" != "1" ]]; then
  EWAS_DATAHUB_LOG="${EWAS_DATAHUB_LOG:-$LOGDIR/ewas_datahub_EWAS_db.log}" \
    bash "$SCRIPT_DIR/post_ewas_datahub_download.sh"
fi
