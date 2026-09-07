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
# Parallelism for missing GSM only (default 8). Override with EWAS_RETRY_JOBS=1 for serial.
JOBS="${EWAS_RETRY_JOBS:-8}"

if [[ ! -f "$MANIFEST" ]]; then
  printf 'Manifest not found: %s\nRun: uv run python scripts/summarize_ewas_db_download_failures.py\n' "$MANIFEST" >&2
  exit 1
fi

log() { printf '%s %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a "$LOG"; }

fetch_one() {
  local study="$1" file="$2"
  local out="$DEST/$study/$file"
  mkdir -p "$DEST/$study"
  if [[ -f "$out" && -s "$out" ]]; then
    return 0
  fi
  local url="${HTTP_ROOT}/${study}/${file}"
  printf '%s wget %s/%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$study" "$file" >>"$LOG"
  wget -c --tries=3 --retry-connrefused --waitretry=10 --timeout=60 --read-timeout=120 \
    -q -O "$out" "$url" >>"$LOG" 2>&1 || {
    printf 'WARN: failed %s/%s\n' "$study" "$file" | tee -a "$LOG" >&2
    # Drop empty failed stubs so a later pass retries.
    [[ -f "$out" && ! -s "$out" ]] && rm -f "$out"
  }
}
export -f fetch_one
export DEST HTTP_ROOT LOG

log "retry manifest=$MANIFEST dest=$DEST jobs=$JOBS"

# Build work list of still-missing rows, then fetch in parallel.
work="$(mktemp)"
trap 'rm -f "$work"' EXIT
tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r study file; do
  [[ -z "$study" || -z "$file" ]] && continue
  out="$DEST/$study/$file"
  if [[ -f "$out" && -s "$out" ]]; then
    continue
  fi
  printf '%s\t%s\n' "$study" "$file"
done >"$work"

n_todo="$(wc -l <"$work" | tr -d ' ')"
log "missing_to_fetch=$n_todo"
if [[ "$n_todo" -gt 0 ]]; then
  # Null-delimited study/file pairs (GSM paths have no spaces).
  while IFS=$'\t' read -r study file; do
    printf '%s\0%s\0' "$study" "$file"
  done <"$work" | xargs -0 -P "$JOBS" -n 2 bash -c 'fetch_one "$1" "$2"' _ || true
fi

log "retry pass complete; run post hook or summarize script to refresh manifest"
if [[ "${EWAS_DATAHUB_SKIP_POST_HOOK:-0}" != "1" ]]; then
  EWAS_DATAHUB_LOG="${EWAS_DATAHUB_LOG:-$LOGDIR/ewas_datahub_EWAS_db.log}" \
    bash "$SCRIPT_DIR/post_ewas_datahub_download.sh"
fi
