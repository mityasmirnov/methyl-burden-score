#!/usr/bin/env bash
# Train study-grouped deepMAT flat baselines on Hub-derived matrices + write reports.
set -euo pipefail
cd "$(dirname "$0")/.."
# shellcheck disable=SC1091
source scripts/activate_data_environment.sh

export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
MAX_LOCI="${MAX_LOCI:-8000}"
LOGDIR="${MBS_ARTIFACT_ROOT}/logs/train"
mkdir -p "$LOGDIR"
TS="$(date -u +%Y%m%dT%H%M%SZ)"

train_one() {
  local cfg="$1"
  local run_id="$2"
  echo "=== train ${run_id} config=${cfg} max_loci=${MAX_LOCI} ==="
  uv run mbs train flat \
    --config "$cfg" \
    --run-id "$run_id" \
    --device cuda \
    --max-loci "$MAX_LOCI" \
    2>&1 | tee "$LOGDIR/${run_id}_${TS}.log"
}

train_one configs/experiment/stage0_flat_hub_age.yaml stage0-hub-age-studyholdout-v1
train_one configs/experiment/stage0_flat_hub_tissue.yaml stage0-hub-tissue-studyholdout-v1
train_one configs/experiment/stage0_flat_hub_blood.yaml stage0-hub-blood-studyholdout-v1
train_one configs/experiment/stage0_flat_hub_brain.yaml stage0-hub-brain-studyholdout-v1
# Disease profile pack still downloading; age-pack overlap lacks case labels.
# train_one configs/experiment/stage0_flat_hub_disease.yaml stage0-hub-disease-studyholdout-v1

echo "=== writing benchmark reports ==="
uv run python scripts/write_hub_real_benchmark_report.py

echo "=== hub real train+report finished ==="
