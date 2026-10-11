#!/usr/bin/env bash
# Run ONE job on GPU0 under the sharing protocol (docs/TODO_PIPELINE.md, "GPU0
# SHARING PROTOCOL"), instead of each agent hand-rolling a queue script.
#
#   scripts/run_gpu0_job.sh <job-name> -- <command ...>
#
# Refuses to start if:
#   * scratch/GPU0_OWNER exists and its pid is alive   (someone already owns it)
#   * another run_gpu0_queue_* script is running       (the 2026-09-28 overlap)
#   * GPU0 has less than MIN_FREE_GIB free             (default 8)
# A lock whose pid is dead is reported and replaced (stale after a crash/kill).
# The lock is removed on exit, including on SIGTERM/SIGINT.
#
# Advisory only: cowrd's gpu_wait does not read the lock. It makes ownership
# visible and lets any agent see who to message. Never touches other users' jobs.
set -uo pipefail
cd "$(dirname "$0")/.."

JOB="${1:-}"; shift || true
[[ "${1:-}" == "--" ]] && shift
if [[ -z "$JOB" || $# -eq 0 ]]; then
  echo "usage: $0 <job-name> -- <command ...>" >&2; exit 2
fi

LOCK=scratch/GPU0_OWNER
MIN_FREE_GIB="${MIN_FREE_GIB:-8}"
mkdir -p scratch/logs

if [[ -f "$LOCK" ]]; then
  old_pid=$(awk -F= '/^pid=/{print $2}' "$LOCK")
  if [[ -n "$old_pid" ]] && kill -0 "$old_pid" 2>/dev/null; then
    echo "REFUSED: GPU0 is owned (pid $old_pid). Lock:" >&2; cat "$LOCK" >&2; exit 3
  fi
  echo "note: replacing STALE lock (pid ${old_pid:-?} is dead):" >&2; cat "$LOCK" >&2
fi

if pgrep -af "run_gpu0_queue" | grep -v pgrep | grep -q .; then
  echo "REFUSED: another run_gpu0_queue_* script is running:" >&2
  pgrep -af "run_gpu0_queue" | grep -v pgrep >&2; exit 3
fi

free_gib=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits -i 0 | awk '{printf "%d", $1/1024}')
if (( free_gib < MIN_FREE_GIB )); then
  echo "REFUSED: GPU0 has only ${free_gib} GiB free (< ${MIN_FREE_GIB})." >&2; exit 3
fi

TS=$(date +%Y%m%d_%H%M%S)
LOG="scratch/logs/gpu0_${JOB}_${TS}.log"
cat > "$LOCK" <<EOF
owner=methyl-burden-score
job=${JOB}
pid=$$
started=$(date -Is)
log=${LOG}
yield=message the methyl-burden-score Claude session or ask the user; job stops within minutes
EOF
CHILD=""
cleanup() { rm -f "$LOCK"; echo "[$(date -Is)] lock released (${JOB})" | tee -a "$LOG"; }
trap cleanup EXIT
# bash defers traps while a FOREGROUND command runs, so a TERM to this wrapper would
# not act until the job finished and "yield within minutes" would be false. Run the
# job in the background and `wait`, which a trapped signal interrupts immediately;
# then forward the signal to the job.
on_signal() {
  echo "[$(date -Is)] signal received, stopping ${JOB} (pid ${CHILD})" | tee -a "$LOG"
  [[ -n "$CHILD" ]] && kill -TERM "$CHILD" 2>/dev/null
  [[ -n "$CHILD" ]] && wait "$CHILD" 2>/dev/null
  exit 143
}
trap on_signal TERM INT

export CUDA_VISIBLE_DEVICES=0
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export PYTHONUNBUFFERED=1
source scripts/activate_data_environment.sh 2>/dev/null || true

echo "[$(date -Is)] START ${JOB}: $*" | tee -a "$LOG"
"$@" >> "$LOG" 2>&1 &
CHILD=$!
wait "$CHILD"
rc=$?
echo "[$(date -Is)] END ${JOB} rc=${rc}" | tee -a "$LOG"
exit $rc
