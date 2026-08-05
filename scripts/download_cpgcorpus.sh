#!/usr/bin/env bash
# Sync full CpGCorpus into $MBS_DATA_ROOT/raw/cpgcorpus.
# Does nothing until explicitly invoked. Logs to $MBS_ARTIFACT_ROOT/logs/downloads.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

if ! command -v aws >/dev/null 2>&1; then
  printf 'ERROR: aws CLI is required for CpGCorpus download.\n' >&2
  exit 1
fi

TARGET="$MBS_DATA_ROOT/raw/cpgcorpus"
LOGDIR="$MBS_ARTIFACT_ROOT/logs/downloads"
mkdir -p "$TARGET" "$LOGDIR" "$MBS_SCRATCH_ROOT/downloads"

STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
LOGFILE="$LOGDIR/cpgcorpus_full_${STAMP}.log"

export AWS_DEFAULT_REGION="${AWS_DEFAULT_REGION:-us-east-1}"
export AWS_REGION="${AWS_REGION:-$AWS_DEFAULT_REGION}"

exec > >(tee -a "$LOGFILE") 2>&1

printf 'Log file: %s\n' "$LOGFILE"
printf 'Syncing CpGCorpus into %s\n' "$TARGET"
printf 'Requester-pays S3; AWS credentials must be configured.\n'

aws s3 sync \
  s3://cpgpt-lucascamillo-public/data/cpgcorpus/raw \
  "$TARGET" \
  --request-payer requester \
  --region "$AWS_DEFAULT_REGION"

printf 'CpGCorpus sync finished: %s\n' "$TARGET"
