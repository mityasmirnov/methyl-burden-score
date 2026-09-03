# GEO backfill pilot summary

- Generated: `2026-09-03T09:22:12Z`
- GEO parquet GSM in: **23217**
- Catalog samples touched (EWAS_db-only): **20621**
- Hub-skipped GSM: **2210**
- GEO GSM not in catalog: **386**
- Phenotype rows added: **32967**
- Samples with ≥1 observed GEO phenotype: **18599**
- Samples with `metadata_json.geo`: **20621**

## Invariants

- GEO phenotype rows on Hub pack-membership GSM: **0** (must be 0)
- Atlas blobs on `sample.metadata_json`: **0** (must be 0)
- `sample_source_membership` includes `geo_metadata_backfill`: **False** (must be false)

## Phenotypes by id

| phenotype_id | rows | observed | unique GSM |
| --- | ---: | ---: | ---: |
| `age` | 7652 | 7652 | 7652 |
| `disease` | 2058 | 2058 | 2058 |
| `sex` | 12729 | 12729 | 12729 |
| `tissue` | 10528 | 10528 | 10528 |

## Label status

| phenotype_id | label_status | n |
| --- | --- | ---: |
| `age` | `observed` | 7652 |
| `disease` | `control` | 2058 |
| `sex` | `observed` | 12729 |
| `tissue` | `observed` | 10528 |

## Eligibility (`source_family=geo_metadata_backfill`)

| phenotype_id | n | cases | controls | core | aux |
| --- | ---: | ---: | ---: | --- | --- |
| `age` | 7652 | None | None | True | True |
| `disease` | 2058 | 0 | 2058 | False | False |
  - `disease` not core: need ≥200 cases, ≥200 controls, ≥3 studies (unknown≠control)
| `sex` | 12729 | None | None | False | True |
  - `sex` not core: sex is auxiliary biological / QC, not a core burden target
| `tissue` | 10528 | None | None | True | True |

## Per study (merge)

| study_id | samples touched | phenotype rows |
| --- | ---: | ---: |
| `GSE105018` | 1230 | 1230 |
| `GSE109379` | 1104 | 1104 |
| `GSE130051` | 1501 | 1501 |
| `GSE140686` | 1504 | 1504 |
| `GSE145361` | 1889 | 4708 |
| `GSE147740` | 1128 | 3384 |
| `GSE157131` | 1218 | 2436 |
| `GSE185920` | 1471 | 2942 |
| `GSE197678` | 2922 | 5844 |
| `GSE210255` | 1394 | 1394 |
| `GSE224124` | 1107 | 1107 |
| `GSE270375` | 994 | 0 |
| `GSE55763` | 1841 | 5523 |
| `GSE56046` | 290 | 290 |
| `GSE68379` | 1028 | 0 |

## Census delta

- Unique GSM before: **150725**
- Unique GSM after: **150947**
- GEO observed phenotype rows before: **32967**
- GEO observed phenotype rows after: **32967**

GEO backfill does not add `sample` rows (EWAS_db scan does). Unique-GSM movement is EWAS_db mirror growth, not GEO. Phenotype-row movement is the GEO delta.

## Operator notes

- GEO rows are omitted for Hub GSM (Hub wins).
- Atlas enrichment stays on study.metadata_json only.
- Disease/cancer rows need explicit case/control tokens.
- Diagnosis-only text stays in metadata_json.
- Training heads read Hub pack Parquet, not geo_metadata_backfill.
