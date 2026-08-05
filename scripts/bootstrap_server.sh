#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${MBS_PROJECT_ROOT:-/data/projects/methyl-burden-score}"
DATA_ROOT="${MBS_DATA_ROOT:-/data/datasets/methyl-burden-score}"
SCRATCH_ROOT="${MBS_SCRATCH_ROOT:-/data/scratch/methyl-burden-score}"
CACHE_ROOT="${MBS_CACHE_ROOT:-/data/cache/methyl-burden-score}"
ARTIFACT_ROOT="${MBS_ARTIFACT_ROOT:-/data/artifacts/methyl-burden-score}"
TOOLS_ROOT="${MBS_TOOLS_ROOT:-/data/tools}"

for path in \
  "$PROJECT_ROOT" \
  "$DATA_ROOT/raw" \
  "$DATA_ROOT/staging" \
  "$DATA_ROOT/canonical/catalog/tables" \
  "$DATA_ROOT/canonical/matrices" \
  "$DATA_ROOT/canonical/annotations" \
  "$DATA_ROOT/canonical/graphs" \
  "$DATA_ROOT/canonical/static_features" \
  "$SCRATCH_ROOT/tmp" \
  "$CACHE_ROOT" \
  /data/cache/uv \
  /data/cache/pip \
  /data/cache/huggingface \
  /data/cache/torch \
  /data/cache/xdg \
  /data/cache/numba \
  /data/cache/matplotlib \
  "$ARTIFACT_ROOT/runs" \
  "$ARTIFACT_ROOT/checkpoints" \
  "$ARTIFACT_ROOT/scores" \
  "$ARTIFACT_ROOT/reports" \
  "$ARTIFACT_ROOT/wandb" \
  "$TOOLS_ROOT"; do
  mkdir -p "$path"
done

if [[ ! -f "$PROJECT_ROOT/.env" ]]; then
  cp "$PROJECT_ROOT/.env.example" "$PROJECT_ROOT/.env"
  printf 'Created %s/.env from template\n' "$PROJECT_ROOT"
fi

# shellcheck disable=SC1091
source "$PROJECT_ROOT/scripts/activate_data_environment.sh"

if ! command -v uv >/dev/null 2>&1; then
  UV_INSTALL_DIR="$TOOLS_ROOT/uv/bin" UV_NO_MODIFY_PATH=1 \
    curl -LsSf https://astral.sh/uv/install.sh | sh
  export PATH="$TOOLS_ROOT/uv/bin:$PATH"
fi

cd "$PROJECT_ROOT"
uv venv --python 3.11 .venv
uv sync --all-groups --all-extras

printf '\nBootstrap completed.\n'
printf 'Project:   %s\n' "$PROJECT_ROOT"
printf 'Data:      %s\n' "$DATA_ROOT"
printf 'Scratch:   %s\n' "$SCRATCH_ROOT"
printf 'Cache:     %s\n' "$CACHE_ROOT"
printf 'Artifacts: %s\n' "$ARTIFACT_ROOT"
printf '\nActivate future shells with:\n  source %s/scripts/activate_data_environment.sh\n' "$PROJECT_ROOT"
