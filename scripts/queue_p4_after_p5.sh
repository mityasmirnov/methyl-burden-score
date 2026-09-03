#!/usr/bin/env bash
# Wait for P5-G-max-explicit 3/3 folds, then launch P4-G-explicit alone on GPU 0.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

LOGDIR="$MBS_ARTIFACT_ROOT/logs/train"
mkdir -p "$LOGDIR"
LOG="$LOGDIR/stage0_7g_queue_p4_after_p5.log"
P5_RUN="stage0-7g-gene-probe-P5-G-max-explicit"
P5_ROOT="$MBS_ARTIFACT_ROOT/runs/$P5_RUN"

log() { printf '[%s] %s\n' "$(date -u +%H:%M:%S)" "$*" | tee -a "$LOG"; }

fold_count() {
  local n=0
  for i in 0 1 2; do
    [[ -f "$P5_ROOT/fold_$i/metrics.json" ]] && n=$((n + 1))
  done
  printf '%s' "$n"
}

log "queue: waiting for $P5_RUN 3/3 folds before P4-G-explicit"
while [[ "$(fold_count)" -lt 3 ]]; do
  log "P5-G-max-explicit: $(fold_count)/3 folds"
  sleep 60
done
log "P5 complete; launching P4-G-explicit on GPU 0 (share=1, 2GiB encoder-parity reserve)"
export CUDA_VISIBLE_DEVICES=0
export DEVICE=cuda
export MBS_CASCADE_GPU_SHARE=1
export MBS_CASCADE_GPU_RESERVED_MIB=2048

CUDA_VISIBLE_DEVICES=0 ARMS="P4-G" DEVICE=cuda \
  MBS_CASCADE_GPU_SHARE=1 MBS_CASCADE_GPU_RESERVED_MIB=2048 \
  bash "$REPO_ROOT/scripts/train_7g_gene_only_probe_background.sh" \
  >>"$LOGDIR/stage0_7g_P4G_after_p5.log" 2>&1
log "P4-G launcher exited $?"
