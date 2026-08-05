#!/usr/bin/env bash
# Download ALL EWAS DataHub public data into $MBS_DATA_ROOT/raw/ewas_datahub:
#   1) Baseline Data packs (HTTP *_v1.zip)
#   2) All Data tree (FTP EWAS_db/)
# Sources: https://ngdc.cncb.ac.cn/ewas/datahub/download
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

TARGET="$MBS_DATA_ROOT/raw/ewas_datahub"
BASELINE_DIR="$TARGET/baseline"
ALLDATA_DIR="$TARGET/EWAS_db"
mkdir -p "$BASELINE_DIR" "$ALLDATA_DIR" \
  "$MBS_ARTIFACT_ROOT/logs/downloads" "$MBS_SCRATCH_ROOT/downloads"

HTTP_BASE="https://download.cncb.ac.cn/ewas/datahub/download"
FTP_ALL="ftp://download.big.ac.cn/ewas/datahub/EWAS_db/"
HTTP_ALL="https://download.cncb.ac.cn/ewas/datahub/EWAS_db/"
FTP_BASELINE="ftp://download.big.ac.cn/ewas/datahub/download/"

BASELINE_FILES=(
  tissue_methylation_v1.zip
  sample_tissue_methylation_v1.zip
  brain_methylation_v1.zip
  sample_brain_methylation_v1.zip
  blood_methylation_v1.zip
  sample_blood_methylation_v1.zip
  sex_methylation_v1.zip
  sample_sex_methylation_v1.zip
  age_methylation_v1.zip
  sample_age_methylation_v1.zip
  ancestry_category_methylation_v1.zip
  sample_ancestry_category_methylation_v1.zip
  bmi_methylation_v1.zip
  sample_bmi_methylation_v1.zip
  cancer_methylation_v1.zip
  sample_cancer_methylation_v1.zip
  disease_methylation_v1.zip
  sample_disease_methylation_v1.zip
)

cat > "$TARGET/SOURCE.txt" <<EOF
source: EWAS DataHub @ EWAS Open Platform
portal: https://ngdc.cncb.ac.cn/ewas/datahub/download
policy: all public data (All Data FTP + Baseline packs)
all_data_ftp: ${FTP_ALL}
all_data_http: ${HTTP_ALL}
baseline_ftp: ${FTP_BASELINE}
baseline_http: ${HTTP_BASE}
gmqn_pmid: 35069703
downloaded: $(date -u +%Y-%m-%dT%H:%M:%SZ)
EOF

printf '=== EWAS DataHub: Baseline Data (HTTP) -> %s ===\n' "$BASELINE_DIR"
for name in "${BASELINE_FILES[@]}"; do
  printf '  %s\n' "$name"
  # Keep legacy flat path if a previous run already wrote there.
  if [[ -f "$TARGET/$name" && ! -f "$BASELINE_DIR/$name" ]]; then
    mv "$TARGET/$name" "$BASELINE_DIR/$name"
  fi
  wget -c -O "$BASELINE_DIR/$name" "$HTTP_BASE/$name"
done

printf '=== EWAS DataHub: All Data (HTTP mirror) -> %s ===\n' "$ALLDATA_DIR"
printf 'Primary: %s\n' "$HTTP_ALL"
printf 'FTP fallback (often stalls here): %s\n' "$FTP_ALL"
# Prefer HTTP; FTP from this host frequently hangs on connect.
set +e
wget --continue --recursive --no-parent --no-host-directories \
  --cut-dirs=3 \
  --reject 'index.html*' \
  --directory-prefix="$ALLDATA_DIR" \
  "$HTTP_ALL"
http_status=$?
if [[ "$http_status" -ne 0 ]]; then
  printf 'WARN: HTTP All Data mirror exited %s; trying FTP...\n' "$http_status" >&2
  wget --continue --recursive --no-parent --no-host-directories \
    --cut-dirs=3 \
    --directory-prefix="$ALLDATA_DIR" \
    "$FTP_ALL"
  ftp_status=$?
  if [[ "$ftp_status" -ne 0 ]]; then
    printf 'WARN: All Data FTP mirror exited %s. Baseline packs are still usable.\n' "$ftp_status" >&2
    printf 'Manual fallback: FileZilla -> %s -> %s\n' "$FTP_ALL" "$ALLDATA_DIR" >&2
  fi
fi
set -e

printf 'EWAS DataHub download finished: %s\n' "$TARGET"
