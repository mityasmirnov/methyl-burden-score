# Hub membership, platform coverage, and DataHub metadata census

**Status:** implemented  
**Parent:** [`data-infrastructure-improvements.md`](data-infrastructure-improvements.md) §2  
**Related:** [`EWAS_DATA.md`](../EWAS_DATA.md), [`DATA_CONTRACT.md`](../DATA_CONTRACT.md),
[`geo-metadata-backfill-ewas-db.md`](geo-metadata-backfill-ewas-db.md)

## Scope and acceptance

Make Hub vs All-Data membership first-class, fill `platform_id` from the official
CNCB repository census (GEO GPL fallback), and ingest curated All-Data metadata
(platform / tissue / disease + sparse field bag) without replacing Hub packs as
training phenotype SoT.

**Done when:**

1. `sample_lane_flags` / `study_lane_flags` (+ views) on refresh; membership report.
2. `make fetch-ewas-datahub-census` writes census Parquet + vocabs; platform totals
   match official 450K/850K/935K (± API).
3. Refresh merges census → platform / tissue / sample_type / `metadata_json.datahub`;
   Hub packs win on phenotype overlap.
4. Unit tests for flags, platform map, merge, cache-only fetch (no network in CI).

## Locked decisions

| Choice | Decision |
|--------|----------|
| Normative All-Data census | CNCB `GET /ewas/datahub/repository/basic?field=platform&val=…` |
| Hub training phenotypes | Unchanged; packs win |
| Platform map | `450K→HM450`, `850K→EPIC`, `935K→EPICv2` |
| Disease rows from census | Only when `sample type=disease tissue`; unknown ≠ control |
| Skip merge | `MBS_SKIP_DATAHUB_CENSUS=1` |

## Commands

```bash
source scripts/activate_data_environment.sh

# Full official census (~180k rows; resumable cache)
make fetch-ewas-datahub-census
# or: uv run python scripts/fetch_ewas_datahub_repository_census.py --from-cache-only

MBS_SKIP_ATLAS_SEED=1 make catalog-refresh-release
```

## Artifacts

| Path | Role |
|------|------|
| `$MBS_CACHE_ROOT/ewas_datahub_repository/{450K,850K,935K}/page_*.json` | API page cache |
| `$MBS_DATA_ROOT/canonical/phenotypes/ewas_datahub_sample_census.parquet` | One row per sample |
| `…/ewas_datahub_sample_census.manifest.json` | Totals + checksum |
| `…/ewas_datahub_field_vocabulary.json` | Union of field keys |
| `…/ewas_datahub_{tissue,disease}_vocab.parquet` | Value counts |
| `catalog/tables/sample_lane_flags.parquet` | Hub vs EWAS_db per sample |
| `catalog/tables/study_lane_flags.parquet` | Aggregate per study |
| `reports/inspection/deepmat_data_v1/hub_vs_ewas_db_membership.{json,md}` | Membership census |
| `reports/inspection/deepmat_data_v1/datahub_metadata_census/summary.{json,md}` | Official vs local |

## Module

[`src/mbs/datahub_census.py`](../../src/mbs/datahub_census.py) — fetch, lane flags, merge, reports.  
SQL: [`sql/014_sample_lane_flags.sql`](../../sql/014_sample_lane_flags.sql).

Refresh writes a **compact** `metadata_json.datahub` bag (platform / tissue /
sample_type / disease / age / sex / project_id / catalog_platform_id). The full
sparse 269-key rows remain in the census Parquet for inspection.

## Non-goals

- Encoder features from DataHub fields.
- Replacing Hub pack phenotype SoT.
- Waiting for an FTP metadata dump (does not exist on download.cncb).

## Observed census discrepancies (2026-09-07 fetch)

| Metric | Official UI | Local census |
|--------|------------:|-------------:|
| Samples | 180 317 | **180 317** (exact) |
| 450K / 850K / 935K | homepage split | **105 085 / 74 187 / 1 045** |
| Tissues/Cells | 1 413 | 1 415 |
| Diseases | 842 | 847 |
| Fields | 296 | 269 (union of keys present on any row; sparse columns absent from all pages) |

Tissue/disease deltas are small (± a few values). Field count is lower because the
repository API returns sparse objects — keys never populated on any downloaded
page never appear in the union vocabulary.
