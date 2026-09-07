#!/usr/bin/env bash
# CPU chain: wait for fold-0 panel, classical-only, then remaining folds.
# Does NOT launch neural Stage B. Safe while GPU 0 is on Milestone 10.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
# shellcheck disable=SC1091
source "$ROOT/scripts/activate_data_environment.sh"
cd "$ROOT"
# Cap BLAS threads — 65k-col SGD thrashing with default all-core OpenBLAS.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-8}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-8}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-8}"
PANEL_DIR="$ROOT/reports/inspection/stage0_7g_prime_matched_probe/fold_panels"
LOG="$ROOT/scratch/logs/stage_b_cpu_chain.log"
mkdir -p "$ROOT/scratch/logs" "$PANEL_DIR"

log() { printf '[%s] %s\n' "$(date -Is)" "$*" | tee -a "$LOG"; }

log "waiting for fold_0_panel.json (panels-only fold 0 job)"
while [[ ! -f "$PANEL_DIR/fold_0_panel.json" ]]; do
  sleep 60
done
log "fold 0 panel present; classical-only fold 0"
uv run python -u scripts/run_7g_prime_stage_b.py --device cpu --classical-only --folds 0 \
  >>"$LOG" 2>&1
log "panels-only folds 1,2"
uv run python -u scripts/run_7g_prime_stage_b.py --device cpu --panels-only --folds 1,2 \
  >>"$LOG" 2>&1
log "classical-only folds 1,2"
uv run python -u scripts/run_7g_prime_stage_b.py --device cpu --classical-only --folds 1,2 \
  >>"$LOG" 2>&1
log "CPU chain complete"
