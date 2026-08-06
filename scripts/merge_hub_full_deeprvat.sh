#!/usr/bin/env bash
# Merge full Hub age/tissue/sex matrices + phenotype table (Milestone 5d).
set -euo pipefail
cd "$(dirname "$0")/.."
# shellcheck disable=SC1091
source scripts/activate_data_environment.sh

OUT_ID="${OUT_ID:-matrix-hub-age-tissue-sex-full-v1}"
if [[ -d "data/canonical/matrices/${OUT_ID}" ]]; then
  echo "skip existing merge matrix_id=${OUT_ID}"
  exit 0
fi

uv run mbs phenotypes build-multitask-table \
  --matrix-id "$OUT_ID" \
  --age-matrix-id matrix-hub-age-full-v1 \
  --tissue-matrix-id matrix-hub-tissue-full-v1 \
  --sex-matrix-id matrix-hub-sex-full-v1 \
  --phenotype-table data/canonical/phenotypes/sample_phenotype_table_age_tissue_sex_full_v1.parquet \
  --tissue-ontology data/canonical/phenotypes/tissue_ontology_age_tissue_sex_full_v1.yaml \
  --sex-ontology data/canonical/phenotypes/sex_ontology_v1.yaml

echo "MERGE_DONE matrix_id=${OUT_ID}"
