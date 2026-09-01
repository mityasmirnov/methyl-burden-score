#!/usr/bin/env bash
# Background 7G′ Stage A gene-only probe: P*-G + classical -G + orphan ablation.
# Survives Cursor/session disconnect via nohup. Re-run to resume (--skip-if-done on cascade).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

CONFIG="${CONFIG:-$REPO_ROOT/configs/experiment/stage0_7g_gene_only_probe.yaml}"
REPORT_DIR="${REPORT_DIR:-$REPO_ROOT/reports/inspection/stage0_7g_gene_only_probe}"
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
DEVICE="${DEVICE:-cuda}"

LOGDIR="$MBS_ARTIFACT_ROOT/logs/train"
mkdir -p "$LOGDIR"
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
NOHUP_LOG="$LOGDIR/stage0_7g_gene_only_${STAMP}.log"
PID_FILE="$LOGDIR/stage0_7g_gene_only_${STAMP}.pid"
LATEST_LOG="$LOGDIR/stage0_7g_gene_only_latest.log"
LATEST_PID="$LOGDIR/stage0_7g_gene_only_latest.pid"

export PYTHONUNBUFFERED=1
nohup uv run python "$REPO_ROOT/scripts/run_7g_gene_only_probe.py" \
  --config "$CONFIG" \
  --device "$DEVICE" \
  >"$NOHUP_LOG" 2>&1 &
PID=$!
printf '%s\n' "$PID" >"$PID_FILE"
ln -sfn "$NOHUP_LOG" "$LATEST_LOG"
ln -sfn "$PID_FILE" "$LATEST_PID"

printf 'Started background 7G′ Stage A gene-only probe\n'
printf '  config: %s\n' "$CONFIG"
printf '  pid:    %s\n' "$PID"
printf '  log:    %s\n' "$NOHUP_LOG"
printf '  pidfile:%s\n' "$PID_FILE"
printf '  GPU:    CUDA_VISIBLE_DEVICES=%s\n' "$CUDA_VISIBLE_DEVICES"
printf 'Status: bash scripts/status_7g_gene_only_probe.sh\n'
printf 'Poll:   tail -f %s\n' "$LATEST_LOG"
