#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${MBS_PROJECT_ROOT:-/data/projects/methyl-burden-score}"
cd "$PROJECT_ROOT"

status=0

check_variable() {
  local name="$1"
  local value="${!name:-}"
  if [[ -z "$value" ]]; then
    printf 'WARN: %s is not set\n' "$name"
    return
  fi
  if [[ "$value" != /data/* ]]; then
    printf 'ERROR: %s=%s is outside /data\n' "$name" "$value" >&2
    status=1
  fi
}

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
  check_variable "$variable"
done

for forbidden in "$HOME/.cache" "$HOME/.local/share/uv" "$HOME/.cache/huggingface"; do
  if [[ -e "$forbidden" ]]; then
    size="$(du -sh "$forbidden" 2>/dev/null | cut -f1 || true)"
    printf 'WARN: existing home artifact %s (%s); inspect before moving or deleting\n' "$forbidden" "${size:-unknown}"
  fi
done

if grep -RInE '(^|[=:[:space:]])(~|/home/[^/]+)/(cache|data|scratch|artifacts|models|checkpoints)' \
  --exclude-dir=.git \
  --exclude-dir=.venv \
  --exclude='*.md' \
  .; then
  printf 'ERROR: found runtime paths under home in tracked implementation files\n' >&2
  status=1
fi

exit "$status"
