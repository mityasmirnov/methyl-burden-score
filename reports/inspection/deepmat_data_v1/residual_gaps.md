# Residual catalog / DataHub / platform gaps

- Generated: `2026-09-28T09:33:07Z`
- Release: `deepmat-data-v1`
- Applied assay platform fill this run: **True** (set 18 / considered 19)

## Catalog (release DuckDB)

- Samples: **185,851**
- Studies: **1,860**
- Assay files: **187,324**

| Platform | Studies |
| --- | ---: |
| `HM450` | 1162 |
| `EPIC` | 636 |
| `(null)` | 42 |
| `EPICv2` | 20 |

## EWAS_db on disk

- Study dirs: **1989** (with `*.txt`: **1860**, empty: **129**)
- Sample files: **187,324** (`GSM*.txt`: **173,166**)

## Catalog vs DataHub census

- Census unique sample IDs (comma-split): **180,634**
- Catalog samples: **185,851**
- Overlap: **169,367**
- In census, not in catalog: **11,267**
- In catalog, not in census: **16,484**

Top studies for census-not-in-catalog:

| Study | Missing census IDs |
| --- | ---: |
| `GSE169156` | 2052 |
| `CPTAC-3` | 585 |
| `GSE152026` | 486 |
| `GSE179847` | 479 |
| `GSE154915` | 470 |
| `GSE167885` | 411 |
| `GSE168739` | 407 |
| `GSE159907` | 316 |
| `GSE166844` | 217 |
| `GSE172365` | 206 |
| `GSE183015` | 174 |
| `GSE184269` | 167 |
| `GSE171140` | 160 |
| `GSE179414` | 157 |
| `GSE164822` | 153 |
| `GSE183798` | 145 |
| `GSE164056` | 143 |
| `GSE156299` | 137 |
| `GSE174818` | 128 |
| `GSE180061` | 121 |

### Why the census gap remains

Empty local dirs with remote_txt=0 cannot be filled from the Hub mirror. Partial census>disk studies are Hub-side incomplete, not local download debt.

- Remote probe empty-on-both-sides: **128**
- Remote recoverable (delta>0): **0**

## Study `platform_id` nulls

- Still null: **42** studies (**9,833** samples)

| Reason | Studies |
| --- | ---: |
| `mixed_census` | 39 |
| `no_single_platform_source` | 2 |
| `single_census_unapplied` | 1 |

- `mixed_census`: Hub census lists ≥2 Illumina platforms for the study — study-level `platform_id` correctly stays null (`metadata_json.datahub.platforms` holds the set).
- `assay_probe_count`: fillable from EWAS_db `*.txt` line counts (HM450/EPIC/EPICv2 bands).
- `no_single_platform_source`: no census/GEO/assay agreement.

## Notes

- Planned empty-dir refill lists (`configs/data/ewas_db_refill_*.txt`) are fully applied; remaining empty dirs have **no** Hub remote `*.txt`.
- Non-GSM namespaces (TCGA/CPTAC/ArrayExpress/…) are mirrored where Hub hosts files.
- Re-run: `uv run python scripts/write_residual_gaps_report.py [--apply-platform-fill]`
