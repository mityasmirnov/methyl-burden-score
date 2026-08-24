#!/usr/bin/env bash
# Cursor Cloud environment helper (idempotent, safe to re-run).
#
# The MBS filesystem policy (AGENTS.md, src/mbs/paths.py) requires MBS_ROOT and
# every data/cache/scratch/artifact root to resolve under /data. On the Cursor
# Cloud VM the git working tree lives at /workspace, so this script bind-mounts
# the working tree at /data/projects/methyl-burden-score. A bind mount (unlike a
# symlink) keeps Path.resolve() under /data, which the path-policy checks and the
# test-suite fixtures rely on.
#
# It does NOT install dependencies (use `uv sync`) and does NOT export shell
# variables (source scripts/activate_data_environment.sh for that). It needs sudo
# only for the one-time `mkdir /data` and the bind mount.
set -euo pipefail

WORKTREE="${MBS_WORKTREE:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
PROJECT_ROOT="${MBS_ROOT:-/data/projects/methyl-burden-score}"

# 1. Ensure /data exists and is owned by the current user (one-time, needs sudo).
if [[ ! -d /data ]]; then
  sudo mkdir -p /data
  sudo chown "$(id -u):$(id -g)" /data
fi
mkdir -p "$(dirname "$PROJECT_ROOT")"
mkdir -p "$PROJECT_ROOT"

# 2. Bind-mount the working tree so the repo genuinely lives under /data.
if ! mountpoint -q "$PROJECT_ROOT"; then
  sudo mount --bind "$WORKTREE" "$PROJECT_ROOT"
fi

# 3. Seed .env in the working tree if absent (git-ignored; project-local paths
#    that place data/cache/scratch/artifacts under the /data project root).
if [[ ! -f "$WORKTREE/.env" ]]; then
  cp "$WORKTREE/.env.example" "$WORKTREE/.env"
fi

echo "cloud_data_setup complete:"
echo "  worktree     = $WORKTREE"
echo "  project_root = $PROJECT_ROOT (bind mount)"
echo "Next: cd $PROJECT_ROOT && source scripts/activate_data_environment.sh && uv sync --all-groups --extra training --extra analysis --frozen"
