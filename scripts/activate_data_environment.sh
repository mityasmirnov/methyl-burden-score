#!/usr/bin/env bash
# Source this file; do not execute it in a subshell.

set -a
if [[ -f /data/projects/methyl-burden-score/.env ]]; then
  # shellcheck disable=SC1091
  source /data/projects/methyl-burden-score/.env
fi
set +a

export MBS_PROJECT_ROOT="${MBS_PROJECT_ROOT:-/data/projects/methyl-burden-score}"
export MBS_DATA_ROOT="${MBS_DATA_ROOT:-/data/datasets/methyl-burden-score}"
export MBS_SCRATCH_ROOT="${MBS_SCRATCH_ROOT:-/data/scratch/methyl-burden-score}"
export MBS_CACHE_ROOT="${MBS_CACHE_ROOT:-/data/cache/methyl-burden-score}"
export MBS_ARTIFACT_ROOT="${MBS_ARTIFACT_ROOT:-/data/artifacts/methyl-burden-score}"
export MBS_DOCKER_ROOT="${MBS_DOCKER_ROOT:-/data/docker}"

export UV_CACHE_DIR="${UV_CACHE_DIR:-/data/cache/uv}"
export PIP_CACHE_DIR="${PIP_CACHE_DIR:-/data/cache/pip}"
export HF_HOME="${HF_HOME:-/data/cache/huggingface}"
export HUGGINGFACE_HUB_CACHE="${HUGGINGFACE_HUB_CACHE:-$HF_HOME/hub}"
export TRANSFORMERS_CACHE="${TRANSFORMERS_CACHE:-$HF_HOME/transformers}"
export TORCH_HOME="${TORCH_HOME:-/data/cache/torch}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-/data/cache/xdg}"
export NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-/data/cache/numba}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/data/cache/matplotlib}"

export TMPDIR="${TMPDIR:-$MBS_SCRATCH_ROOT/tmp}"
export TEMP="$TMPDIR"
export TMP="$TMPDIR"
export WANDB_DIR="${WANDB_DIR:-$MBS_ARTIFACT_ROOT/wandb}"
export WANDB_CACHE_DIR="${WANDB_CACHE_DIR:-$MBS_CACHE_ROOT/wandb}"

export PATH="/data/tools/uv/bin:$MBS_PROJECT_ROOT/.venv/bin:$PATH"
export PYTHONPATH="$MBS_PROJECT_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONUNBUFFERED="${PYTHONUNBUFFERED:-1}"
export TOKENIZERS_PARALLELISM="${TOKENIZERS_PARALLELISM:-false}"

mkdir -p \
  "$MBS_DATA_ROOT" \
  "$MBS_SCRATCH_ROOT/tmp" \
  "$MBS_CACHE_ROOT" \
  "$MBS_ARTIFACT_ROOT" \
  "$UV_CACHE_DIR" \
  "$PIP_CACHE_DIR" \
  "$HF_HOME" \
  "$TORCH_HOME" \
  "$XDG_CACHE_HOME" \
  "$NUMBA_CACHE_DIR" \
  "$MPLCONFIGDIR" \
  "$WANDB_DIR" \
  "$WANDB_CACHE_DIR"

for variable in \
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
