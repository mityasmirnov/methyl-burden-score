#!/usr/bin/env bash
# Download EWAS Atlas batch association exports into $MBS_DATA_ROOT/raw/ewas_atlas
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

TARGET="$MBS_DATA_ROOT/raw/ewas_atlas"
mkdir -p "$TARGET" "$MBS_ARTIFACT_ROOT/logs/downloads" "$MBS_SCRATCH_ROOT/downloads"

BASE="https://download.cncb.ac.cn/ewas"
FILES=(
  EWAS_Atlas_associations.tsv
  EWAS_Atlas_studies.tsv
  EWAS_Atlas_cohorts.tsv
  EWAS_Atlas_probe_annotations.tsv
  EWAS_trait_trait_logP.txt
)

printf 'Downloading EWAS Atlas batch files into %s\n' "$TARGET"
for name in "${FILES[@]}"; do
  printf '  %s\n' "$name"
  wget -c -O "$TARGET/$name" "$BASE/$name"
done

printf 'EWAS Atlas download finished: %s\n' "$TARGET"
