#!/usr/bin/env bash
# Background GEO crawl chain: species rebuild → fetch next 100 → catalog refresh
# → eligibility + geo-dev. Does not mutate ATS or launch training.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1091
source "$ROOT/scripts/activate_data_environment.sh"
cd "$ROOT"
mkdir -p scratch/logs
LOG="scratch/logs/geo_crawl_chain_$(date -u +%Y%m%dT%H%M%SZ).log"
exec > >(tee -a "$LOG") 2>&1
log() { printf '%s %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*"; }

log "=== GEO crawl chain start ==="
log "log=$LOG"

log "1/5 rebuild species on cached batch parquet"
uv run python scripts/fetch_geo_sample_metadata.py \
  --from-cache-only \
  --studies-file configs/data/geo_backfill_cached_batch_gse.txt

log "2/5 fetch next 100 GSE (SOFT download + parquet merge)"
uv run python scripts/fetch_geo_sample_metadata.py \
  --studies-file configs/data/geo_backfill_next_gse.txt

log "3/5 catalog refresh (Hub∩GEO fill-missing + GEO merge)"
MBS_SKIP_ATLAS_SEED=1 make catalog-refresh-release

log "4/5 eligibility-by-study audit"
uv run python scripts/write_geo_eligibility_by_study.py

log "5/5 rebuild deepmat-data-geo-dev-v1 phenotype arms"
uv run python scripts/build_geo_dev_release.py

log "=== GEO crawl chain done ==="
log "ATS not mutated; training not launched"
printf '%s\n' "$LOG" > scratch/logs/geo_crawl_chain_latest.path
