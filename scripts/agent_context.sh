#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${MBS_PROJECT_ROOT:-/data/projects/methyl-burden-score}"
cd "$PROJECT_ROOT"

printf '=== methyl-burden-score agent context ===\n'
printf 'project_root: %s\n' "$PROJECT_ROOT"
printf 'git_branch:   %s\n' "$(git branch --show-current)"
printf 'git_commit:   %s\n' "$(git rev-parse HEAD)"
printf 'python:       %s\n' "$(python --version 2>&1 || true)"
printf 'uv:           %s\n' "$(uv --version 2>&1 || true)"
printf '\n=== working tree ===\n'
git status --short

printf '\n=== top-level project map ===\n'
find . \
  -maxdepth 3 \
  -not -path './.git*' \
  -not -path './.venv*' \
  -not -path './vendor/*/.git*' \
  -not -path './data*' \
  -not -path './artifacts*' \
  -not -path './scratch*' \
  | sort \
  | sed -n '1,240p'

printf '\n=== reference repositories ===\n'
if [[ -f .gitmodules ]]; then
  git submodule status || true
else
  printf 'not installed; run make references\n'
fi

printf '\n=== authoritative documents ===\n'
printf '%s\n' \
  AGENTS.md \
  docs/ARCHITECTURE.md \
  docs/DATA_CONTRACT.md \
  docs/ANNOTATION_GRAPH.md \
  docs/EXPERIMENT_PROTOCOL.md \
  docs/WORKSPACE.md

printf '\n=== data roots ===\n'
for variable in MBS_DATA_ROOT MBS_SCRATCH_ROOT MBS_CACHE_ROOT MBS_ARTIFACT_ROOT; do
  value="${!variable:-not-set}"
  printf '%-20s %s\n' "$variable" "$value"
done
