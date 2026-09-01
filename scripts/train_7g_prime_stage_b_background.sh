#!/usr/bin/env bash
# Background 7G′ Stage B matched-panel benchmark.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

CONFIG="${CONFIG:-$REPO_ROOT/configs/experiment/stage0_7g_prime_stage_b.yaml}"
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
DEVICE="${DEVICE:-cuda}"

LOGDIR="$MBS_ARTIFACT_ROOT/logs/train"
mkdir -p "$LOGDIR"
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
NOHUP_LOG="$LOGDIR/stage0_7g_prime_stage_b_${STAMP}.log"
PID_FILE="$LOGDIR/stage0_7g_prime_stage_b_${STAMP}.pid"

export PYTHONUNBUFFERED=1
nohup uv run python "$REPO_ROOT/scripts/run_7g_prime_stage_b.py" \
  --config "$CONFIG" \
  --device "$DEVICE" \
  >"$NOHUP_LOG" 2>&1 &
PID=$!
printf '%s\n' "$PID" >"$PID_FILE"
printf 'Started background 7G′ Stage B (pid=%s log=%s)\n' "$PID" "$NOHUP_LOG"
