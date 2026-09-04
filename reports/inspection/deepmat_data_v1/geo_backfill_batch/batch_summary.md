# GEO batch-50 summary

- Generated: `2026-09-04T15:17:44Z`
- Parquet GSM: **45762** across **50** GSE
- Catalog samples touched: **38020**
- Hub-skipped: **7077**
- Phenotype rows added this merge: **74971**
- Catalog GEO phenotype rows: **74971**
- Hub overlap GEO rows: **0** (must be 0)
- Multi-study GSM: **0**

## Tissue map (touched GSM)

- {'ambiguous': 0, 'empty': 16899, 'mapped': 11670, 'unmapped': 9451}

## Phenotypes by id (merge)

- {'age': 19233, 'cancer': 1470, 'disease': 5216, 'sex': 27931, 'tissue': 21121}

## Label status

| phenotype_id | label_status | n |
| --- | --- | ---: |
| `age` | `observed` | 19233 |
| `cancer` | `case` | 1239 |
| `cancer` | `control` | 231 |
| `disease` | `case` | 1456 |
| `disease` | `control` | 3760 |
| `sex` | `observed` | 27931 |
| `tissue` | `observed` | 21121 |

## Eligibility (GEO family)

| phenotype_id | n | cases | controls | core |
| --- | ---: | ---: | ---: | --- |
| `age` | 19233 | None | None | True |
| `cancer` | 1470 | 1239 | 231 | True |
| `disease` | 5216 | 1456 | 3760 | True |
| `sex` | 27931 | None | None | False |
  - not core: sex is auxiliary biological / QC, not a core burden target
| `tissue` | 21121 | None | None | True |

## Notes

- Training heads still read Hub pack Parquet, not geo_metadata_backfill.
- Do not train GEO disease/cancer until eligibility clears cases+controls.
- Do not mutate frozen ATS; GEO-enriched training release is separate/not built.
