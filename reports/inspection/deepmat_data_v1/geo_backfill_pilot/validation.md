# GEO repaired-pilot validation (15 GSE)

- Generated: `2026-09-04T14:35:06Z`
- Parquet GSM: **23217**
- Multi-study GSM (membership): **0**
- Fetch conflict samples: **0**

## Global tissue map

| status | n |
| --- | ---: |
| `empty` | 11434 |
| `mapped` | 7289 |
| `unmapped` | 4494 |

## Global age units (after conversion to years)

| age_unit | n |
| --- | ---: |
| `years` | 9435 |

- Age years min/median/max: **5.99** / **41.02121834** / **83.0** (n=9435)

## Catalog (after merge)

- GEO phenotype rows: **33961**
- Disease/cancer **cases**: **994** (must not train disease head if 0)
- GEO rows on Hub membership GSM: **0** (must be 0)

## Tissue coverage by study

| study_id | gsm | mapped | unmapped | empty | top unmapped |
| --- | ---: | ---: | ---: | ---: | --- |
| `GSE105018` | 1658 | 0 | 0 | 1658 | — |
| `GSE109379` | 1104 | 0 | 1104 | 0 | brain tumor (1104) |
| `GSE130051` | 1501 | 0 | 1501 | 0 | Medulloblastoma (1501) |
| `GSE140686` | 1505 | 0 | 1505 | 0 | sarcoma (1505) |
| `GSE145361` | 1889 | 1889 | 0 | 0 | — |
| `GSE147740` | 1129 | 0 | 0 | 1129 | — |
| `GSE157131` | 1218 | 1218 | 0 | 0 | — |
| `GSE185920` | 1471 | 1471 | 0 | 0 | — |
| `GSE197678` | 2922 | 0 | 0 | 2922 | — |
| `GSE210255` | 1394 | 0 | 0 | 1394 | — |
| `GSE224124` | 1107 | 0 | 0 | 1107 | — |
| `GSE270375` | 1378 | 0 | 384 | 994 | Meningioma (384) |
| `GSE55763` | 2711 | 2711 | 0 | 0 | — |
| `GSE56046` | 1202 | 0 | 0 | 1202 | — |
| `GSE68379` | 1028 | 0 | 0 | 1028 | — |

## Operator notes

- Training heads read Hub pack Parquet, not geo_metadata_backfill.
- GEO improves catalog / eligibility only until a separate geo-dev release.
- Disease/cancer cases must be >0 with eligibility cutoffs before any disease head.
- Do not mutate frozen ATS (matrix-hub-age-tissue-sex-full-v1).
