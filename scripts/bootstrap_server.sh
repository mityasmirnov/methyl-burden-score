#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

if [[ -f "$REPO_ROOT/.env" ]]; then
  # shellcheck disable=SC1091
  source "$REPO_ROOT/.env"
fi

export MBS_ROOT="${MBS_ROOT:-${MBS_PROJECT_ROOT:-$REPO_ROOT}}"
export MBS_PROJECT_ROOT="${MBS_PROJECT_ROOT:-$MBS_ROOT}"
export MBS_DATA_ROOT="${MBS_DATA_ROOT:-$MBS_ROOT/data}"
export MBS_SCRATCH_ROOT="${MBS_SCRATCH_ROOT:-$MBS_ROOT/scratch}"
export MBS_CACHE_ROOT="${MBS_CACHE_ROOT:-$MBS_ROOT/cache}"
export MBS_ARTIFACT_ROOT="${MBS_ARTIFACT_ROOT:-$MBS_ROOT/artifacts}"
export MBS_DOCKER_ROOT="${MBS_DOCKER_ROOT:-$MBS_ROOT/docker}"

mkdir -p \
  "$MBS_DATA_ROOT/raw/cpgcorpus" \
  "$MBS_DATA_ROOT/raw/ewas_datahub" \
  "$MBS_DATA_ROOT/raw/ewas_atlas" \
  "$MBS_DATA_ROOT/raw/manifests" \
  "$MBS_DATA_ROOT/raw/gencode" \
  "$MBS_DATA_ROOT/raw/annotations" \
  "$MBS_DATA_ROOT/staging" \
  "$MBS_DATA_ROOT/staging/infinium_export" \
  "$MBS_DATA_ROOT/canonical/catalog/tables" \
  "$MBS_DATA_ROOT/canonical/matrices" \
  "$MBS_DATA_ROOT/canonical/annotations" \
  "$MBS_DATA_ROOT/canonical/graphs" \
  "$MBS_DATA_ROOT/canonical/static_features" \
  "$MBS_SCRATCH_ROOT/tmp" \
  "$MBS_SCRATCH_ROOT/downloads" \
  "$MBS_CACHE_ROOT/xdg" \
  "$MBS_CACHE_ROOT/huggingface" \
  "$MBS_CACHE_ROOT/torch" \
  "$MBS_CACHE_ROOT/uv" \
  "$MBS_CACHE_ROOT/pip" \
  "$MBS_CACHE_ROOT/numba" \
  "$MBS_CACHE_ROOT/matplotlib" \
  "$MBS_ARTIFACT_ROOT/runs" \
  "$MBS_ARTIFACT_ROOT/checkpoints" \
  "$MBS_ARTIFACT_ROOT/scores" \
  "$MBS_ARTIFACT_ROOT/reports" \
  "$MBS_ARTIFACT_ROOT/wandb" \
  "$MBS_ARTIFACT_ROOT/logs/downloads" \
  "$MBS_DOCKER_ROOT" \
  "$MBS_ROOT/.tools/uv/bin"

if [[ ! -f "$MBS_ROOT/.env" ]]; then
  cp "$MBS_ROOT/.env.example" "$MBS_ROOT/.env"
  printf 'Created %s/.env from template\n' "$MBS_ROOT"
fi

# shellcheck disable=SC1091
source "$MBS_ROOT/scripts/activate_data_environment.sh"

UV_BIN_DIR="$MBS_ROOT/.tools/uv/bin"
mkdir -p "$UV_BIN_DIR"

ensure_project_uv() {
  if [[ -x "$UV_BIN_DIR/uv" ]]; then
    return 0
  fi
  if command -v uv >/dev/null 2>&1; then
    # Relocate a pre-existing install (e.g. $HOME/.local/bin) under /data.
    cp "$(command -v uv)" "$UV_BIN_DIR/uv"
    chmod +x "$UV_BIN_DIR/uv"
    if command -v uvx >/dev/null 2>&1; then
      cp "$(command -v uvx)" "$UV_BIN_DIR/uvx"
      chmod +x "$UV_BIN_DIR/uvx"
    fi
    return 0
  fi
  # Official installer; force project-local destination, then verify.
  UV_INSTALL_DIR="$UV_BIN_DIR" UV_NO_MODIFY_PATH=1 \
    curl -LsSf https://astral.sh/uv/install.sh | sh
  if [[ ! -x "$UV_BIN_DIR/uv" ]] && command -v uv >/dev/null 2>&1; then
    cp "$(command -v uv)" "$UV_BIN_DIR/uv"
    chmod +x "$UV_BIN_DIR/uv"
  fi
  if [[ ! -x "$UV_BIN_DIR/uv" ]]; then
    printf 'ERROR: uv was not installed under %s\n' "$UV_BIN_DIR" >&2
    exit 1
  fi
}

ensure_project_uv
export PATH="$UV_BIN_DIR:$PATH"

cd "$MBS_ROOT"
uv venv --python 3.11 .venv
uv sync --all-groups --all-extras

echo "Bootstrap complete."
echo "MBS_ROOT=$MBS_ROOT"
echo "MBS_DATA_ROOT=$MBS_DATA_ROOT"
echo "MBS_SCRATCH_ROOT=$MBS_SCRATCH_ROOT"
echo "MBS_CACHE_ROOT=$MBS_CACHE_ROOT"
echo "MBS_ARTIFACT_ROOT=$MBS_ARTIFACT_ROOT"
echo "uv=$UV_BIN_DIR/uv"
echo
echo "Activate future shells with:"
echo "  source $MBS_ROOT/scripts/activate_data_environment.sh"
