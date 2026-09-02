#!/usr/bin/env bash
# Regenerate Stage A analysis.md + lock from existing per_arm/*.json (no retraining).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/scripts/activate_data_environment.sh"

REPORT="${REPORT_DIR:-$REPO_ROOT/reports/inspection/stage0_7g_gene_only_probe}"
uv run python "$REPO_ROOT/scripts/write_7g_gene_only_probe_report.py" --report-dir "$REPORT"
printf 'Regenerated %s/analysis.md and lock_recommendation.json\n' "$REPORT"
