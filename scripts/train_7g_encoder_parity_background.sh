#!/usr/bin/env bash
# Background 7G′ optional encoder parity (Flat-G + Hier-G on gene_cols).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

CONFIG="${CONFIG:-$REPO_ROOT/configs/experiment/stage0_7g_encoder_parity.yaml}"
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
DEVICE="${DEVICE:-cuda}"
FORCE="${FORCE:-}"

LOGDIR="$MBS_ARTIFACT_ROOT/logs/train"
mkdir -p "$LOGDIR"
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
NOHUP_LOG="$LOGDIR/stage0_7g_encoder_parity_${STAMP}.log"
PID_FILE="$LOGDIR/stage0_7g_encoder_parity_${STAMP}.pid"
LATEST_LOG="$LOGDIR/stage0_7g_encoder_parity_latest.log"
LATEST_PID="$LOGDIR/stage0_7g_encoder_parity_latest.pid"

LOCK="$REPO_ROOT/reports/inspection/stage0_7g_gene_only_probe/lock_recommendation.json"
if [[ -f "$LOCK" ]] && ! grep -q '"recommend_encoder_parity": true' "$LOCK"; then
  if [[ "$FORCE" != "1" ]]; then
    printf 'Skip encoder parity: lock recommends against (set FORCE=1 to override)\n'
    exit 0
  fi
fi

EXTRA=()
if [[ "$FORCE" == "1" ]]; then
  EXTRA+=(--force)
fi

export PYTHONUNBUFFERED=1
nohup uv run python "$REPO_ROOT/scripts/run_7g_encoder_parity_probe.py" \
  --config "$CONFIG" \
  --device "$DEVICE" \
  --skip-if-done \
  "${EXTRA[@]}" \
  >"$NOHUP_LOG" 2>&1 &
PID=$!
printf '%s\n' "$PID" >"$PID_FILE"
ln -sfn "$NOHUP_LOG" "$LATEST_LOG"
ln -sfn "$PID_FILE" "$LATEST_PID"

printf 'Started background 7G′ encoder parity (pid=%s log=%s GPU=%s)\n' "$PID" "$NOHUP_LOG" "$CUDA_VISIBLE_DEVICES"
