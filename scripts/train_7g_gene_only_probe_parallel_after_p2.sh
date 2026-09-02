#!/usr/bin/env bash
# Wait for P2-G-explicit 3/3 folds, stop sequential runner, launch P4 + P5-max in parallel.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

LOGDIR="$MBS_ARTIFACT_ROOT/logs/train"
LOG="$LOGDIR/stage0_7g_parallel_orchestrator.log"
P2_RUN="stage0-7g-gene-probe-P2-G-explicit"
P2_ROOT="$MBS_ARTIFACT_ROOT/runs/$P2_RUN"
SEQUENTIAL_PID="${SEQUENTIAL_PID:-258787}"

log() { printf '[%s] %s\n' "$(date -u +%H:%M:%S)" "$*" | tee -a "$LOG"; }

fold_count() {
  local n=0
  for i in 0 1 2; do
    [[ -f "$P2_ROOT/fold_$i/metrics.json" ]] && n=$((n + 1))
  done
  printf '%s' "$n"
}

log "orchestrator started; waiting for $P2_RUN 3/3 folds"
BATCH_MARKER="$REPO_ROOT/artifacts/logs/train/stage0_7g_batched_cascade.ready"
while [[ "$(fold_count)" -lt 3 ]]; do
  log "P2-G-explicit: $(fold_count)/3 folds"
  sleep 30
done
while [[ ! -f "$BATCH_MARKER" ]]; do
  log "P2 complete; waiting for batched cascade marker at $BATCH_MARKER"
  sleep 15
done
log "batched cascade ready"
log "P2-G-explicit complete; stopping sequential runner pid=$SEQUENTIAL_PID"

if kill -0 "$SEQUENTIAL_PID" 2>/dev/null; then
  kill "$SEQUENTIAL_PID" 2>/dev/null || true
  sleep 5
  kill -9 "$SEQUENTIAL_PID" 2>/dev/null || true
fi
# Stop any cascade child still tied to the old sequential arm queue (not parallel jobs).
pkill -f "run_7g_gene_only_probe.py.*--arm P2-G --arm P4-G" 2>/dev/null || true
sleep 2

log "launching P4-G and P5-G-max on GPU 0 in parallel (gpu_share=2, 2GiB encoder parity reserve)"
export CUDA_VISIBLE_DEVICES=0
export DEVICE=cuda
export MBS_CASCADE_GPU_SHARE=2
export MBS_CASCADE_GPU_RESERVED_MIB=2048

CUDA_VISIBLE_DEVICES=0 ARMS="P4-G" DEVICE=cuda MBS_CASCADE_GPU_SHARE=2 MBS_CASCADE_GPU_RESERVED_MIB=2048 \
  nohup bash "$REPO_ROOT/scripts/train_7g_gene_only_probe_background.sh" \
  >>"$LOGDIR/stage0_7g_P4G_parallel.log" 2>&1 &
P4_PID=$!
log "P4-G pid=$P4_PID log=$LOGDIR/stage0_7g_P4G_parallel.log"

CUDA_VISIBLE_DEVICES=0 ARMS="P5-G-max" DEVICE=cuda MBS_CASCADE_GPU_SHARE=2 MBS_CASCADE_GPU_RESERVED_MIB=2048 \
  nohup bash "$REPO_ROOT/scripts/train_7g_gene_only_probe_background.sh" \
  >>"$LOGDIR/stage0_7g_P5max_parallel.log" 2>&1 &
P5_PID=$!
log "P5-G-max pid=$P5_PID log=$LOGDIR/stage0_7g_P5max_parallel.log"

wait "$P4_PID" && log "P4-G finished" || log "P4-G exited non-zero"
wait "$P5_PID" && log "P5-G-max finished" || log "P5-G-max exited non-zero"

log "optional P5-G-mean (gated)"
CUDA_VISIBLE_DEVICES=0 ARMS="P5-G-mean" DEVICE=cuda \
  bash "$REPO_ROOT/scripts/train_7g_gene_only_probe_background.sh" \
  >>"$LOGDIR/stage0_7g_P5mean_gated.log" 2>&1 || log "P5-G-mean skipped or failed"

log "regenerating report"
uv run python "$REPO_ROOT/scripts/write_7g_gene_only_probe_report.py" >>"$LOG" 2>&1
log "orchestrator done"
