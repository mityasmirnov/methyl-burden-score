#!/usr/bin/env bash
# Convert full Hub packs for Milestone 5d (age/tissue/sex) and 7B (remaining six).
# Skip existing matrix dirs. Age/tissue/sex freezes must not be overwritten.
set -euo pipefail
cd "$(dirname "$0")/.."
# shellcheck disable=SC1091
source scripts/activate_data_environment.sh
LOGDIR="${MBS_ARTIFACT_ROOT}/logs/matrix_convert"
mkdir -p "$LOGDIR"
TS="$(date -u +%Y%m%dT%H%M%SZ)"
LOG="$LOGDIR/hub_full_packs_${TS}.log"
exec > >(tee -a "$LOG") 2>&1

echo "=== Hub full-pack converts starting ${TS} ==="

export_family() {
  local family="$1"
  local parquet="data/canonical/phenotypes/${family}_sample_info.parquet"
  if [[ -f "$parquet" ]]; then
    echo "--- sample-info already present: ${parquet}"
    return 0
  fi
  echo "--- export sample-info family=${family}"
  make export-ewas-sample-info FAMILY="$family"
}

run_full() {
  local family="$1"
  local matrix_id="$2"
  if [[ -f "data/canonical/matrices/${matrix_id}/conversion_stats.json" ]]; then
    echo "--- skip existing matrix_id=${matrix_id}"
    return 0
  fi
  echo "--- convert family=${family} matrix_id=${matrix_id} --all-studies"
  uv run mbs matrix convert-pack \
    --phenotype-family "$family" \
    --matrix-id "$matrix_id" \
    --all-studies \
    --platform-id HM450
  echo "CONVERT_FAMILY_DONE family=${family} matrix_id=${matrix_id}"
}

# 5d freezes (skip if present)
export_family age
export_family tissue
export_family sex
run_full age "matrix-hub-age-full-v1"
run_full tissue "matrix-hub-tissue-full-v1"
run_full sex "matrix-hub-sex-full-v1"

# 7B remaining packs (smaller first so a restart yields evidence sooner)
for family in ancestry bmi brain blood cancer disease; do
  export_family "$family"
  run_full "$family" "matrix-hub-${family}-full-v1"
done

echo "--- build virtual hub pack index + overlap check"
uv run mbs matrix index-hub-packs \
  --check-overlap \
  --report-dir reports/inspection/stage0_7b_hub_matrices

echo "CONVERT_DONE"
echo "=== Hub full-pack converts finished ==="
ls -la data/canonical/matrices/matrix-hub-*-full-v1/ || true
