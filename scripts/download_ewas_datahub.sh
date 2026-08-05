#!/usr/bin/env bash
# Download EWAS DataHub reference methylation packs into $MBS_DATA_ROOT/raw/ewas_datahub
# Sources: https://ngdc.cncb.ac.cn/ewas/datahub/download
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

TARGET="$MBS_DATA_ROOT/raw/ewas_datahub"
mkdir -p "$TARGET" "$MBS_ARTIFACT_ROOT/logs/downloads" "$MBS_SCRATCH_ROOT/downloads"

BASE="https://download.cncb.ac.cn/ewas/datahub/download"
FILES=(
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

printf 'Downloading EWAS DataHub packs into %s\n' "$TARGET"
printf 'These archives are large (tens of GB); resume is enabled via wget -c.\n'
for name in "${FILES[@]}"; do
  printf '  %s\n' "$name"
  wget -c -O "$TARGET/$name" "$BASE/$name"
done

printf 'EWAS DataHub download finished: %s\n' "$TARGET"
