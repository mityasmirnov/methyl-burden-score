#!/usr/bin/env bash
# Launch a CpGCorpus download in the background via nohup.
# Does not start any download unless this script is invoked explicitly.
#
# Usage:
#   scripts/download_cpgcorpus_background.sh gse [extra args...]
#   scripts/download_cpgcorpus_background.sh full
#
# Examples:
#   scripts/download_cpgcorpus_background.sh gse
#   scripts/download_cpgcorpus_background.sh gse GSE116992 GSE35069
#   scripts/download_cpgcorpus_background.sh full
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

MODE="${1:-}"
if [[ -z "$MODE" ]]; then
  printf 'Usage: %s <gse|full> [args...]\n' "$(basename "$0")" >&2
  exit 1
fi
shift || true

if ! command -v aws >/dev/null 2>&1; then
  printf 'ERROR: aws CLI is required for CpGCorpus download.\n' >&2
  exit 1
fi

LOGDIR="$MBS_ARTIFACT_ROOT/logs/downloads"
mkdir -p "$LOGDIR" "$MBS_DATA_ROOT/raw/cpgcorpus" "$MBS_SCRATCH_ROOT/downloads"

STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
case "$MODE" in
  gse|cpgcorpus_gse|cpgcorpus-gse)
    WORKER="$REPO_ROOT/scripts/download_cpgcorpus_gse.sh"
    NOHUP_LOG="$LOGDIR/cpgcorpus_gse_nohup_${STAMP}.log"
    PID_FILE="$LOGDIR/cpgcorpus_gse_nohup_${STAMP}.pid"
    ;;
  full|cpgcorpus)
    WORKER="$REPO_ROOT/scripts/download_cpgcorpus.sh"
    NOHUP_LOG="$LOGDIR/cpgcorpus_full_nohup_${STAMP}.log"
    PID_FILE="$LOGDIR/cpgcorpus_full_nohup_${STAMP}.pid"
    ;;
  *)
    printf 'ERROR: unknown mode %s (expected gse or full)\n' "$MODE" >&2
    exit 1
    ;;
esac

nohup bash "$WORKER" "$@" >"$NOHUP_LOG" 2>&1 &
PID=$!
printf '%s\n' "$PID" >"$PID_FILE"

printf 'Started background CpGCorpus download\n'
printf '  mode: %s\n' "$MODE"
printf '  pid:  %s\n' "$PID"
printf '  log:  %s\n' "$NOHUP_LOG"
printf '  pidfile: %s\n' "$PID_FILE"
printf 'Monitor with: tail -f %s\n' "$NOHUP_LOG"
