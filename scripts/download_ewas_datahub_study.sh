#!/usr/bin/env bash
# Mirror one or more EWAS DataHub All Data studies from
#   https://download.cncb.ac.cn/ewas/datahub/EWAS_db/{STUDY}/
# into $MBS_DATA_ROOT/raw/ewas_datahub/EWAS_db/{STUDY}/
#
# Usage:
#   bash scripts/download_ewas_datahub_study.sh GSE35069
#   bash scripts/download_ewas_datahub_study.sh GSE35069 GSE125367
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

if [[ "$#" -lt 1 ]]; then
  printf 'Usage: %s STUDY [STUDY ...]\n' "$(basename "$0")" >&2
  exit 1
fi

TARGET="$MBS_DATA_ROOT/raw/ewas_datahub"
LOGDIR="$MBS_ARTIFACT_ROOT/logs/downloads"
mkdir -p "$TARGET/EWAS_db" "$LOGDIR" "$MBS_SCRATCH_ROOT/downloads"

HTTP_ROOT="https://download.cncb.ac.cn/ewas/datahub"

list_hrefs() {
  local url="$1"
  curl -fsSL --max-time 120 "$url" \
    | grep -oE 'href="[^"]+"' \
    | sed 's/^href="//;s/"$//' \
    | grep -Ev '^\.\./|^/|^javascript:|^#|^https?:|[()]'
}

download_study() {
  local study="$1"
  local root_url="${HTTP_ROOT}/EWAS_db/${study}/"
  local dest="$TARGET/EWAS_db/${study}"
  mkdir -p "$dest"
  printf '=== EWAS_db/%s -> %s ===\n' "$study" "$dest"
  mapfile -t files < <(list_hrefs "$root_url" | grep -Ev '/$' | grep -E '^GSM[0-9]+\.txt$' || true)
  if [[ "${#files[@]}" -eq 0 ]]; then
    printf 'ERROR: no files listed for %s at %s\n' "$study" "$root_url" >&2
    return 1
  fi
  printf 'Found %s files\n' "${#files[@]}"
  local i=0
  for f in "${files[@]}"; do
    i=$((i + 1))
    printf '  [%s/%s] %s\n' "$i" "${#files[@]}" "$f"
    wget -c --tries=3 --retry-connrefused --waitretry=10 --timeout=60 --read-timeout=120 -q \
      -O "$dest/$f" "${root_url}${f}" || {
      printf 'WARN: failed %s/%s\n' "$study" "$f" >&2
    }
  done
  local n_local
  n_local="$(find "$dest" -type f -name '*.txt' | wc -l | tr -d ' ')"
  printf 'Downloaded %s txt files for %s\n' "$n_local" "$study"
}

LOGFILE="$LOGDIR/ewas_datahub_study_$(date -u +%Y%m%dT%H%M%SZ).log"
{
  printf 'start: %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  for study in "$@"; do
    download_study "$study"
  done
  printf 'finished: %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
} 2>&1 | tee -a "$LOGFILE"

printf 'EWAS DataHub study download finished (log: %s)\n' "$LOGFILE"
