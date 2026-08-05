#!/usr/bin/env bash
# Download EPICv2 reannotated manifest (Zenodo) into $MBS_DATA_ROOT/raw/manifests
# Record: https://doi.org/10.5281/zenodo.14933468 (latest file version 20704849)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

TARGET="$MBS_DATA_ROOT/raw/manifests/epicv2"
mkdir -p "$TARGET" "$MBS_ARTIFACT_ROOT/logs/downloads" "$MBS_SCRATCH_ROOT/downloads"

# Recommended Stage 0 file from Zenodo record 20704849
ZENODO_API="https://zenodo.org/api/records/20704849/files"
FILES=(
  EPICv2_reannotated_manifest_v3.0.csv.gz
)

printf 'Downloading EPICv2 reannotated manifest into %s\n' "$TARGET"
for name in "${FILES[@]}"; do
  printf '  %s\n' "$name"
  wget -c -O "$TARGET/$name" "$ZENODO_API/$name/content"
done

cat > "$TARGET/SOURCE.txt" <<EOF
source: Zenodo EPICv2 Re-annotated Manifest v3.0
concept_doi: https://doi.org/10.5281/zenodo.14933468
record: https://zenodo.org/records/20704849
vendor_code: vendor/epicv2_manifest
downloaded: $(date -u +%Y-%m-%dT%H:%M:%SZ)
EOF

printf 'Manifest download finished: %s\n' "$TARGET"
