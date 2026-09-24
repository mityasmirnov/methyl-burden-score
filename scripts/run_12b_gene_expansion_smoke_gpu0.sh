#!/usr/bin/env bash
# Milestone 12b go/no-go smoke on GPU0 (idle while GPU2 runs the Milestone 12
# N-light OOF): full gene-linked column universe (no 65,536-prefix cap) vs.
# the Milestone-10 65k-prefix baseline, for both N-light and P2-G cascade.
# Single fold (fold 0 of hub-nine-pack-3fold-v1), 16 epochs each. Sequential
# on one GPU -- N-light first (cheaper), then cascade.
set -uo pipefail
cd "$(dirname "$0")/.."
source scripts/activate_data_environment.sh

export CUDA_VISIBLE_DEVICES=0
export PYTHONUNBUFFERED=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

log() { echo "[$(date -Is)] $*"; }

REPORT=reports/inspection/stage0_12b_gene_expansion_smoke
mkdir -p "$REPORT"

log "=== 12b gene-expansion smoke: N-light (GPU0) ==="
uv run python -u scripts/run_12b_gene_expansion_nlight_smoke.py
log "=== N-light gene-expansion smoke exit=$? ==="

log "=== 12b gene-expansion smoke: P2-G cascade (GPU0) ==="
uv run mbs train cascade \
  --config configs/experiment/stage0_12b_gene_expansion_cascade_smoke.yaml \
  --run-id stage0-12b-gene-expansion-cascade-smoke-f0 \
  --device cuda \
  --max-folds 1 \
  --skip-if-done \
  --report-dir "$REPORT"
log "=== cascade gene-expansion smoke exit=$? ==="

log "=== 12b gene-expansion smoke done ==="
