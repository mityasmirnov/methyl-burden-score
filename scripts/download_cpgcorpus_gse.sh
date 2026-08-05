#!/usr/bin/env bash
# Sync selected GSE prefixes from CpGCorpus (requester-pays S3).
# Does nothing until explicitly invoked. Logs to $MBS_ARTIFACT_ROOT/logs/downloads.
# Usage:
#   scripts/download_cpgcorpus_gse.sh
#   scripts/download_cpgcorpus_gse.sh GSE42861 GSE87571
#   scripts/download_cpgcorpus_gse.sh --list configs/data/stage0_cpgcorpus_gse_list.txt
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

if ! command -v aws >/dev/null 2>&1; then
  printf 'ERROR: aws CLI is required.\n' >&2
  exit 1
fi

export AWS_DEFAULT_REGION="${AWS_DEFAULT_REGION:-us-east-1}"
export AWS_REGION="${AWS_REGION:-$AWS_DEFAULT_REGION}"

TARGET="$MBS_DATA_ROOT/raw/cpgcorpus"
LOGDIR="$MBS_ARTIFACT_ROOT/logs/downloads"
mkdir -p "$TARGET" "$LOGDIR" "$MBS_SCRATCH_ROOT/downloads"

STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
LOGFILE="$LOGDIR/cpgcorpus_gse_${STAMP}.log"

BUCKET="s3://cpgpt-lucascamillo-public/data/cpgcorpus/raw"
DEFAULT_LIST="$REPO_ROOT/configs/data/stage0_cpgcorpus_gse_list.txt"

GSES=()
if [[ "${1:-}" == "--list" ]]; then
  LIST_FILE="${2:-$DEFAULT_LIST}"
  while IFS= read -r line || [[ -n "$line" ]]; do
    [[ -z "$line" || "$line" =~ ^[[:space:]]*# ]] && continue
    GSES+=("$line")
  done < "$LIST_FILE"
elif [[ "$#" -gt 0 ]]; then
  GSES=("$@")
else
  while IFS= read -r line || [[ -n "$line" ]]; do
    [[ -z "$line" || "$line" =~ ^[[:space:]]*# ]] && continue
    GSES+=("$line")
  done < "$DEFAULT_LIST"
fi

if [[ "${#GSES[@]}" -eq 0 ]]; then
  printf 'ERROR: no GSE IDs provided.\n' >&2
  exit 1
fi

exec > >(tee -a "$LOGFILE") 2>&1

printf 'Log file: %s\n' "$LOGFILE"
printf 'Selective CpGCorpus sync into %s (region=%s)\n' "$TARGET" "$AWS_DEFAULT_REGION"

available=0
missing=0
for gse in "${GSES[@]}"; do
  prefix="${BUCKET}/${gse}/"
  if ! aws s3 ls "$prefix" --request-payer requester --region "$AWS_DEFAULT_REGION" | grep -q .; then
    printf 'MISSING on S3: %s\n' "$gse"
    missing=$((missing + 1))
    continue
  fi
  printf 'Syncing %s\n' "$gse"
  aws s3 sync \
    "$prefix" \
    "$TARGET/$gse/" \
    --request-payer requester \
    --region "$AWS_DEFAULT_REGION"
  available=$((available + 1))
done

printf 'Finished. synced=%s missing=%s target=%s\n' "$available" "$missing" "$TARGET"
if [[ "$available" -eq 0 ]]; then
  exit 2
fi
