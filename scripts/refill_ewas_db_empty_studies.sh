#!/usr/bin/env bash
# Refill empty EWAS_db study dirs from a newline study list (GSM + non-GSM txt).
# Usage:
#   bash scripts/refill_ewas_db_empty_studies.sh configs/data/ewas_db_refill_empty_gse.txt
#   EWAS_REFILL_JOBS=4 bash scripts/refill_ewas_db_empty_studies.sh configs/data/ewas_db_refill_nongsm.txt
#
# Parallelism: EWAS_REFILL_JOBS studies × EWAS_REFILL_FILE_JOBS files/study.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

LIST="${1:?study list required}"
JOBS="${EWAS_REFILL_JOBS:-4}"
FILE_JOBS="${EWAS_REFILL_FILE_JOBS:-4}"
DEST="$MBS_DATA_ROOT/raw/ewas_datahub/EWAS_db"
HTTP_ROOT="https://download.cncb.ac.cn/ewas/datahub/EWAS_db"
LOGDIR="$MBS_ARTIFACT_ROOT/logs/downloads"
mkdir -p "$LOGDIR" "$DEST"
LOG="${EWAS_REFILL_LOG:-$LOGDIR/ewas_db_refill_$(date -u +%Y%m%dT%H%M%SZ).log}"
# Incomplete leftovers from interrupted wget -c (typical complete GSM ≈ 7–14 MiB).
# Typical complete GSM β txt ≈ 7–14 MiB; undersized leftovers from kills re-fetch.
MIN_OK_BYTES="${EWAS_REFILL_MIN_OK_BYTES:-1000000}"

log() { printf '%s %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a "$LOG"; }

list_hrefs() {
  local url="$1"
  curl -fsSL --max-time 120 "$url" \
    | grep -oE 'href="[^"]+"' \
    | sed 's/^href="//;s/"$//' \
    | grep -Ev '^\.\./|^/|^javascript:|^#|^https?:|[()]'
}

download_one_file() {
  local study="$1"
  local f="$2"
  local dest="$DEST/$study"
  local root_url="${HTTP_ROOT}/${study}/"
  local out="$dest/$f"
  local partial="$out.partial"
  if [[ -f "$out" && -s "$out" ]]; then
    local sz
    sz=$(stat -c%s "$out" 2>/dev/null || echo 0)
    if (( sz >= MIN_OK_BYTES )); then
      echo skip
      return 0
    fi
    rm -f "$out"
  fi
  rm -f "$partial"
  if wget -c --tries=3 --retry-connrefused --waitretry=10 --timeout=60 --read-timeout=120 -q \
    -O "$partial" "${root_url}${f}" >>"$LOG" 2>&1; then
    mv -f "$partial" "$out"
    echo ok
  else
    rm -f "$partial"
    printf '%s WARN failed %s/%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$study" "$f" >>"$LOG"
    echo fail
  fi
}

refill_one() {
  local study="$1"
  local dest="$DEST/$study"
  local root_url="${HTTP_ROOT}/${study}/"
  mkdir -p "$dest"
  mapfile -t all_files < <(list_hrefs "$root_url" | grep -Ev '/$' | grep -E '\.txt$' | grep -Ev '^\(\.\+' || true)
  mapfile -t files < <(printf '%s\n' "${all_files[@]:-}" | grep -E '^GSM[0-9]+\.txt$' || true)
  if [[ ${#files[@]} -eq 0 || -z "${files[0]:-}" ]]; then
    mapfile -t files < <(printf '%s\n' "${all_files[@]:-}" | grep -Ev '^\s*$' || true)
  fi
  if [[ ${#files[@]} -eq 0 || -z "${files[0]:-}" ]]; then
    printf '%s SKIP %s (no remote txt)\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$study" >>"$LOG"
    return 0
  fi
  local results
  results=$(
    printf '%s\n' "${files[@]}" \
      | xargs -P "$FILE_JOBS" -I{} bash -c 'download_one_file "$@"' _ "$study" {}
  )
  local ok fail skip
  ok=$(printf '%s\n' "$results" | grep -c '^ok$' || true)
  fail=$(printf '%s\n' "$results" | grep -c '^fail$' || true)
  skip=$(printf '%s\n' "$results" | grep -c '^skip$' || true)
  printf '%s DONE %s ok=%s fail=%s skip=%s listed=%s\n' \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$study" "$ok" "$fail" "$skip" "${#files[@]}" >>"$LOG"
}
export -f refill_one list_hrefs download_one_file log
export DEST HTTP_ROOT LOG FILE_JOBS MIN_OK_BYTES

mapfile -t studies < <(grep -Ev '^\s*(#|$)' "$LIST")
log "refill list=$LIST n_studies=${#studies[@]} jobs=$JOBS file_jobs=$FILE_JOBS dest=$DEST min_ok=$MIN_OK_BYTES"
printf '%s\n' "${studies[@]}" | xargs -P "$JOBS" -I{} bash -c 'refill_one "$@"' _ {}
log "refill pass complete for $LIST"
