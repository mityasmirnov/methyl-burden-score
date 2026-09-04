#!/usr/bin/env bash
# Launch the 7G′ weekend supervisor under nohup with a single flock owner.
# Safe to re-run: supervise_7g_prime_weekend.py refuses a second lock holder.
set -euo pipefail
cd "$(dirname "$0")/.."
source scripts/activate_data_environment.sh

mkdir -p scratch/logs scratch/locks
LOG="scratch/logs/7g_prime_weekend_supervisor.nohup.log"
PIDF="scratch/logs/7g_prime_weekend_launcher.pid"

export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
export PYTHONUNBUFFERED=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

if [[ "${1:-}" == "--status" || "${1:-}" == "--dry-run" ]]; then
  exec uv run python -u scripts/supervise_7g_prime_weekend.py "$@"
fi

# Foreground attach (for debugging)
if [[ "${1:-}" == "--fg" ]]; then
  shift
  exec uv run python -u scripts/supervise_7g_prime_weekend.py "$@"
fi

nohup uv run python -u scripts/supervise_7g_prime_weekend.py "$@" >>"$LOG" 2>&1 &
echo $! >"$PIDF"
echo "started supervisor pid=$(cat "$PIDF") log=$LOG"
echo "status: uv run python scripts/supervise_7g_prime_weekend.py --status"
