#!/usr/bin/env bash
# Source this file; do not execute it in a subshell.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

set -a
if [[ -f "$REPO_ROOT/.env" ]]; then
  # shellcheck disable=SC1091
  source "$REPO_ROOT/.env"
fi
set +a

export MBS_ROOT="${MBS_ROOT:-${MBS_PROJECT_ROOT:-$REPO_ROOT}}"
export MBS_PROJECT_ROOT="${MBS_PROJECT_ROOT:-$MBS_ROOT}"
export MBS_DATA_ROOT="${MBS_DATA_ROOT:-$MBS_ROOT/data}"
export MBS_SCRATCH_ROOT="${MBS_SCRATCH_ROOT:-$MBS_ROOT/scratch}"
export MBS_CACHE_ROOT="${MBS_CACHE_ROOT:-$MBS_ROOT/cache}"
export MBS_ARTIFACT_ROOT="${MBS_ARTIFACT_ROOT:-$MBS_ROOT/artifacts}"
export MBS_DOCKER_ROOT="${MBS_DOCKER_ROOT:-$MBS_ROOT/docker}"

export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$MBS_CACHE_ROOT/xdg}"
export HF_HOME="${HF_HOME:-$MBS_CACHE_ROOT/huggingface}"
export HUGGINGFACE_HUB_CACHE="${HUGGINGFACE_HUB_CACHE:-$HF_HOME/hub}"
export TRANSFORMERS_CACHE="${TRANSFORMERS_CACHE:-$HF_HOME/transformers}"
export TORCH_HOME="${TORCH_HOME:-$MBS_CACHE_ROOT/torch}"
export UV_CACHE_DIR="${UV_CACHE_DIR:-$MBS_CACHE_ROOT/uv}"
export PIP_CACHE_DIR="${PIP_CACHE_DIR:-$MBS_CACHE_ROOT/pip}"
export NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-$MBS_CACHE_ROOT/numba}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-$MBS_CACHE_ROOT/matplotlib}"

export TMPDIR="${TMPDIR:-$MBS_SCRATCH_ROOT/tmp}"
export TEMP="$TMPDIR"
export TMP="$TMPDIR"
export WANDB_DIR="${WANDB_DIR:-$MBS_ARTIFACT_ROOT/wandb}"
export WANDB_CACHE_DIR="${WANDB_CACHE_DIR:-$MBS_CACHE_ROOT/wandb}"

export PATH="$MBS_ROOT/.tools/uv/bin:$MBS_ROOT/.venv/bin:$PATH"
export PYTHONPATH="$MBS_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONUNBUFFERED="${PYTHONUNBUFFERED:-1}"
export TOKENIZERS_PARALLELISM="${TOKENIZERS_PARALLELISM:-false}"

mkdir -p \
  "$MBS_DATA_ROOT" \
  "$MBS_SCRATCH_ROOT/tmp" \
  "$MBS_CACHE_ROOT" \
  "$MBS_ARTIFACT_ROOT" \
  "$MBS_DOCKER_ROOT" \
  "$XDG_CACHE_HOME" \
  "$HF_HOME" \
  "$TORCH_HOME" \
  "$UV_CACHE_DIR" \
  "$PIP_CACHE_DIR" \
  "$NUMBA_CACHE_DIR" \
  "$MPLCONFIGDIR" \
  "$WANDB_DIR" \
  "$WANDB_CACHE_DIR"

for variable in \
  MBS_ROOT \
  MBS_PROJECT_ROOT \
  MBS_DATA_ROOT \
  MBS_SCRATCH_ROOT \
  MBS_CACHE_ROOT \
  MBS_ARTIFACT_ROOT \
  MBS_DOCKER_ROOT \
  UV_CACHE_DIR \
  PIP_CACHE_DIR \
  HF_HOME \
  TORCH_HOME \
  XDG_CACHE_HOME \
  TMPDIR; do
  value="${!variable}"
  if [[ "$value" != /data/* ]]; then
    printf 'ERROR: %s must be under /data, found %s\n' "$variable" "$value" >&2
    return 1 2>/dev/null || exit 1
  fi
done

echo "Environment activated."
echo "MBS_ROOT=$MBS_ROOT"
echo "MBS_DATA_ROOT=$MBS_DATA_ROOT"
echo "MBS_SCRATCH_ROOT=$MBS_SCRATCH_ROOT"
echo "MBS_CACHE_ROOT=$MBS_CACHE_ROOT"
echo "MBS_ARTIFACT_ROOT=$MBS_ARTIFACT_ROOT"
