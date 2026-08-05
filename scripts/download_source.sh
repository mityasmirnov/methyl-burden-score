#!/usr/bin/env bash
# Download a named raw data source under $MBS_DATA_ROOT/raw/<source>.
# Usage: scripts/download_source.sh <source>
# Sources: cpgcorpus | ewas_datahub | ewas_atlas | manifests
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

SOURCE="${1:-}"
if [[ -z "$SOURCE" ]]; then
  printf 'Usage: %s <cpgcorpus|ewas_datahub|ewas_atlas|manifests>\n' "$(basename "$0")" >&2
  exit 1
fi

mkdir -p "$MBS_DATA_ROOT/raw" "$MBS_ARTIFACT_ROOT/logs/downloads" "$MBS_SCRATCH_ROOT/downloads"

case "$SOURCE" in
  cpgcorpus)
    exec bash "$REPO_ROOT/scripts/download_cpgcorpus.sh"
    ;;
  ewas_datahub)
    exec bash "$REPO_ROOT/scripts/download_ewas_datahub.sh"
    ;;
  ewas_atlas)
    exec bash "$REPO_ROOT/scripts/download_ewas_atlas.sh"
    ;;
  manifests)
    exec bash "$REPO_ROOT/scripts/download_manifests.sh"
    ;;
  *)
    printf 'ERROR: unknown source %s\n' "$SOURCE" >&2
    exit 1
    ;;
esac
