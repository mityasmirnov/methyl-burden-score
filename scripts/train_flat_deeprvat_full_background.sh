#!/usr/bin/env bash
# Background train: full Hub DeepRVAT-style age/tissue/sex (Milestone 5d).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

RUN_ID="${RUN_ID:-stage0-flat-deeprvat-age-tissue-sex-full-v1}"
CONFIG="${CONFIG:-$REPO_ROOT/configs/experiment/stage0_flat_deeprvat_full.yaml}"
# Prefer a free GPU; override with CUDA_VISIBLE_DEVICES.
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-1}"

LOGDIR="$MBS_ARTIFACT_ROOT/logs/train"
mkdir -p "$LOGDIR"
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
NOHUP_LOG="$LOGDIR/deeprvat_full_${STAMP}.log"
PID_FILE="$LOGDIR/deeprvat_full_${STAMP}.pid"

nohup uv run mbs train flat \
  --config "$CONFIG" \
  --run-id "$RUN_ID" \
  --device cuda \
  >"$NOHUP_LOG" 2>&1 &
PID=$!
printf '%s\n' "$PID" >"$PID_FILE"

printf 'Started background DeepRVAT full train\n'
printf '  run_id: %s\n' "$RUN_ID"
printf '  pid:    %s\n' "$PID"
printf '  log:    %s\n' "$NOHUP_LOG"
printf '  pidfile:%s\n' "$PID_FILE"
printf '  GPU:    CUDA_VISIBLE_DEVICES=%s\n' "$CUDA_VISIBLE_DEVICES"
printf 'Poll: tail -f %s\n' "$NOHUP_LOG"
printf 'TB:   uv run tensorboard --logdir "$MBS_ARTIFACT_ROOT/runs/%s/tb" --bind_all --port 6007\n' "$RUN_ID"
