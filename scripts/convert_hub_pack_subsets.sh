#!/usr/bin/env bash
# Convert Hub profile pack study-subsets to canonical matrices (real data gate).
set -euo pipefail
cd "$(dirname "$0")/.."
# shellcheck disable=SC1091
source scripts/activate_data_environment.sh
LOGDIR="${MBS_ARTIFACT_ROOT}/logs/matrix_convert"
mkdir -p "$LOGDIR"
TS="$(date -u +%Y%m%dT%H%M%SZ)"
LOG="$LOGDIR/hub_pack_subsets_${TS}.log"
exec > >(tee -a "$LOG") 2>&1

echo "=== Hub pack subset converts starting ${TS} ==="

run_one() {
  local family="$1"
  local studies="$2"
  local matrix_id="$3"
  local max_n="$4"
  echo "--- convert family=${family} matrix_id=${matrix_id} studies=${studies} max_per_study=${max_n}"
  uv run mbs matrix convert-pack \
    --phenotype-family "$family" \
    --study-ids "$studies" \
    --matrix-id "$matrix_id" \
    --max-per-study "$max_n" \
    --platform-id HM450
}

# Age / tissue holdouts for multitask scale-up (max_per_study=100)
run_one age "GSE51032,GSE56105,GSE55763,GSE78874" "matrix-hub-age-studyholdout-v2" 100
run_one tissue "GSE58885,GSE52401,GSE97628,GSE78874,GSE75248" "matrix-hub-tissue-studyholdout-v2" 100

# Keep v1 packs for single-task benchmarks unless already present
if [[ ! -d data/canonical/matrices/matrix-hub-age-studyholdout-v1 ]]; then
  run_one age "GSE51032,GSE56105,GSE55763,GSE78874" "matrix-hub-age-studyholdout-v1" 35
fi
if [[ ! -d data/canonical/matrices/matrix-hub-tissue-studyholdout-v1 ]]; then
  run_one tissue "GSE58885,GSE52401,GSE97628,GSE78874,GSE75248" "matrix-hub-tissue-studyholdout-v1" 35
fi

# Blood pack; train uses tissue column (not cell-fraction phenotype_value)
run_one blood "GSE56105,GSE56046,GSE51032,GSE56581,GSE97628" "matrix-hub-blood-studyholdout-v1" 35

# Brain region multiclass
run_one brain "GSE64509,GSE59457,GSE80970,GSE98203,GSE66351" "matrix-hub-brain-studyholdout-v1" 35

# Disease labels: betas from age pack (disease profile zip still downloading);
# pick age-pack studies that also appear in disease_sample_info.
run_one age "GSE74193,GSE43414,GSE59685" "matrix-hub-disease-from-agepack-v1" 40

echo "=== Rewriting disease matrix sidecar labels from disease_sample_info ==="
uv run python <<'PY'
from pathlib import Path
import pandas as pd

matrix = Path("data/canonical/matrices/matrix-hub-disease-from-agepack-v1")
side = pd.read_parquet(matrix / "sample_phenotypes.parquet")
disease = pd.read_parquet("data/canonical/phenotypes/disease_sample_info.parquet")
by_id = disease.set_index("sample_id")["phenotype_value"].to_dict()
side = side.copy()
side["phenotype_value"] = side["sample_id"].map(by_id)
# keep age numeric if present for diagnostics
side["label_source"] = "disease_sample_info.parquet"
side["beta_source_pack"] = "age_methylation_v1.zip"
missing = int(side["phenotype_value"].isna().sum())
if missing:
    raise SystemExit(f"disease label join missing {missing} samples")
side.to_parquet(matrix / "sample_phenotypes.parquet", index=False)
print("disease sidecar rows", len(side), "label counts:")
print(side["phenotype_value"].fillna("<empty>").astype(str).value_counts().head(20).to_string())
PY

echo "=== Hub pack subset converts finished ==="
ls -la data/canonical/matrices/matrix-hub-*-v*/
