# Phenotype census (deepmat-data-v1)

- Generated: `2026-08-24T12:49:22Z`
- Unique GSM (`sample`): **116113**
- Pack membership row sum: **47843**
- Pack row sum counts Hub membership only; unique GSM also includes EWAS_db-only samples. Pack row sum can exceed unique Hub GSMs when samples appear in multiple packs.

## Pack membership

| Family | Rows | Unique GSM |
| --- | ---: | ---: |
| `age` | 8374 | 8374 |
| `ancestry` | 1380 | 1380 |
| `blood` | 3402 | 3402 |
| `bmi` | 2070 | 2070 |
| `brain` | 1997 | 1997 |
| `cancer` | 10101 | 10101 |
| `disease` | 12218 | 12218 |
| `sex` | 2978 | 2978 |
| `tissue` | 5323 | 5323 |

## Overlap by number of packs

| n_families | n_samples |
| ---: | ---: |
| 1 | 24825 |
| 2 | 5729 |
| 3 | 3177 |
| 4 | 486 |
| 5 | 17 |

## Label conflicts (head)

Rows: 0 (capped at 50 in report).

## EWAS_db ingest

- Local studies: `883` / advertised `1989` (mirror_complete=False)
- Local GSM files: `87153`

Re-run `mbs catalog refresh-release` after more `EWAS_db` study dirs download.

