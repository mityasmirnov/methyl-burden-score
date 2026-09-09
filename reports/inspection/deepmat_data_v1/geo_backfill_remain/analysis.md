# GEO backfill remain — catalog GSE completion

**Date:** 2026-09-09

## Status

GEO sample metadata now covers **all catalog GSE samples** (via primary
`study_id` or multi-series `study_ids` on shared GSM).

| Metric | Value |
|--------|------:|
| Fetch list | 79 catalog GSE previously missing as primary study |
| Download/parse | **79/79** ok, **0** failures (`--view brief`) |
| GEO parquet GSM | **170 338** |
| GEO parquet primary studies | **1 707** |
| Catalog GSE | **1 718** |
| Primary-study gap | **11** (all 100% GSM overlap under related series) |

The 11 “missing” primary IDs (`GSE51057`, `GSE137898`, …) are related/sub-series
of studies already in parquet (e.g. `GSE51057` GSM live under `GSE51032` with
`study_ids` containing both). Sample phenotypes are present; series title/summary
come from `geo_series_metadata` / `study.metadata_json.geo`.

## Phenotype yield (parquet)

| Field | Count |
|-------|------:|
| age non-null | 64 727 |
| sex non-null | 110 430 |
| tissue mapped | 70 244 |
| tissue unmapped | 100 086 |

## EWAS_db assay mirror (separate from GEO labels)

| Field | Value |
|-------|------:|
| Advertised studies | 1 989 |
| Local studies | ~1 695 |
| Local GSM | ~170 641 |
| `mirror_complete` | **false** |

GEO **metadata** backfill for catalog GSE is complete. Residual EWAS_db work is
assay-file mirror gaps / retry failures, not SOFT phenotype crawl.
