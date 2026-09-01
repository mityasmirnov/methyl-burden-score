#!/usr/bin/env bash
# Idempotent 7G cascade tissue probe: P0–P3 + report.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

CONFIG="$REPO_ROOT/configs/experiment/stage0_7g_cascade_tissue_probe.yaml"
DEVICE="cuda"
SKIP_P2=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2 ;;
    --device) DEVICE="$2"; shift 2 ;;
    --skip-p2) SKIP_P2="--skip-p2"; shift ;;
    *) echo "unknown arg: $1" >&2; exit 2 ;;
  esac
done

export PYTHONUNBUFFERED=1
uv run python "$REPO_ROOT/scripts/run_7g_cascade_tissue_probe.py" \
  --config "$CONFIG" \
  --device "$DEVICE" \
  $SKIP_P2
