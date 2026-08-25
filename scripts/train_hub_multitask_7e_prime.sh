#!/usr/bin/env bash
# Milestone 7E′: build virtual Hub cohort + train parameter-matched flat multitask.
# Does not overwrite ATS freeze or v0.1.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
# shellcheck disable=SC1091
source scripts/activate_data_environment.sh

GPU="${CUDA_VISIBLE_DEVICES:-0}"
export CUDA_VISIBLE_DEVICES="$GPU"

echo "==> build virtual Hub multi-store"
uv run mbs matrix build-hub-virtual

echo "==> join Hub union phenotype table"
uv run mbs phenotypes build-hub-union-table

echo "==> train flat Hub multitask (matched encoder)"
uv run mbs train flat \
  --config configs/experiment/stage0_flat_hub_multitask_v1.yaml \
  --run-id stage0-flat-hub-multitask-v1

echo "==> done: artifacts/runs/stage0-flat-hub-multitask-v1/"
