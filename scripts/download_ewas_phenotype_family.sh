#!/usr/bin/env bash
# Download one EWAS Data Hub phenotype family (profile zip + sample-info zip).
# Usage: bash scripts/download_ewas_phenotype_family.sh FAMILY
# FAMILY: age|tissue|disease|cancer|blood|brain|sex|ancestry|bmi
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

FAMILY="${1:-}"
if [[ -z "$FAMILY" ]]; then
  printf 'Usage: %s FAMILY\n' "$0" >&2
  printf 'FAMILY: age tissue disease cancer blood brain sex ancestry bmi\n' >&2
  exit 2
fi

HTTP_ROOT="https://download.cncb.ac.cn/ewas/datahub"
TARGET="$MBS_DATA_ROOT/raw/ewas_datahub/download"
LOGDIR="$MBS_ARTIFACT_ROOT/logs/downloads"
mkdir -p "$TARGET" "$LOGDIR"

case "$FAMILY" in
  age)
    FILES=(age_methylation_v1.zip sample_age_methylation_v1.zip)
    ;;
  tissue)
    FILES=(tissue_methylation_v1.zip sample_tissue_methylation_v1.zip)
    ;;
  disease)
    FILES=(disease_methylation_v1.zip sample_disease_methylation_v1.zip)
    ;;
  cancer)
    FILES=(cancer_methylation_v1.zip sample_cancer_methylation_v1.zip)
    ;;
  blood)
    FILES=(blood_methylation_v1.zip sample_blood_methylation_v1.zip)
    ;;
  brain)
    FILES=(brain_methylation_v1.zip sample_brain_methylation_v1.zip)
    ;;
  sex)
    FILES=(sex_methylation_v1.zip sample_sex_methylation_v1.zip)
    ;;
  ancestry)
    FILES=(ancestry_category_methylation_v1.zip sample_ancestry_category_methylation_v1.zip)
    ;;
  bmi)
    FILES=(bmi_methylation_v1.zip sample_bmi_methylation_v1.zip)
    ;;
  *)
    printf 'Unknown family: %s\n' "$FAMILY" >&2
    exit 2
    ;;
esac

LOG="$LOGDIR/ewas_family_${FAMILY}_$(date -u +%Y%m%dT%H%M%SZ).log"
{
  printf '=== EWAS Data Hub family=%s -> %s ===\n' "$FAMILY" "$TARGET"
  for name in "${FILES[@]}"; do
    printf '  %s\n' "$name"
    if [[ -f "$MBS_DATA_ROOT/raw/ewas_datahub/$name" && ! -f "$TARGET/$name" ]]; then
      mv "$MBS_DATA_ROOT/raw/ewas_datahub/$name" "$TARGET/$name"
    fi
    wget -c -O "$TARGET/$name" "${HTTP_ROOT}/download/${name}"
  done
  printf 'Done family=%s\n' "$FAMILY"
} 2>&1 | tee -a "$LOG"

printf 'Log: %s\n' "$LOG"
