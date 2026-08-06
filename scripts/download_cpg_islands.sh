#!/usr/bin/env bash
# Download UCSC cpIslandExt (hg38) into $MBS_DATA_ROOT/raw/annotations
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

TARGET="$MBS_DATA_ROOT/raw/annotations"
mkdir -p "$TARGET" "$MBS_ARTIFACT_ROOT/logs/downloads" "$MBS_SCRATCH_ROOT/downloads"

URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz"
NAME="cpgIslandExt.hg38.txt.gz"

printf 'Downloading UCSC CpG islands (hg38) into %s\n' "$TARGET"
wget -c -O "$TARGET/$NAME" "$URL"

cat > "$TARGET/SOURCE.txt" <<EOF
source: UCSC Genome Browser goldenPath hg38 cpgIslandExt
file: $NAME
uri: $URL
genome_build: GRCh38
note: table dump; columns include bin, chrom, chromStart, chromEnd (0-based half-open)
shore_bp: 2000
shelf_bp: 4000
downloaded: $(date -u +%Y-%m-%dT%H:%M:%SZ)
EOF

printf 'CpG island download finished: %s\n' "$TARGET/$NAME"
