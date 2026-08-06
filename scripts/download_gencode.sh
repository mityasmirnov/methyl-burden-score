#!/usr/bin/env bash
# Download GENCODE release 38 annotation GTF into $MBS_DATA_ROOT/raw/gencode
# Source: https://www.gencodegenes.org/human/release_38.html
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

TARGET="$MBS_DATA_ROOT/raw/gencode"
mkdir -p "$TARGET" "$MBS_ARTIFACT_ROOT/logs/downloads" "$MBS_SCRATCH_ROOT/downloads"

URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz"
NAME="gencode.v38.annotation.gtf.gz"

printf 'Downloading GENCODE v38 annotation into %s\n' "$TARGET"
wget -c -O "$TARGET/$NAME" "$URL"

cat > "$TARGET/SOURCE.txt" <<EOF
source: GENCODE human release 38
file: $NAME
uri: $URL
genome_build: GRCh38
downloaded: $(date -u +%Y-%m-%dT%H:%M:%SZ)
EOF

printf 'GENCODE download finished: %s\n' "$TARGET/$NAME"
