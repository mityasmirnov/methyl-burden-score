# Data population inventory

- Generated: `2026-09-28T09:25:58Z`
- Release: `deepmat-data-v1`

## Catalog scale

| Samples | Studies | Phenotype rows | Assay files |
| ---: | ---: | ---: | ---: |
| **185,851** | **1,860** | **702,406** | **187,324** |

## Assays / platforms

Platform is study-level (`study.platform_id`); counts are samples under those studies.

| Platform | Samples |
| --- | ---: |
| `HM450` | 113,658 |
| `EPIC` | 59,801 |
| `(null)` | 11,271 |
| `EPICv2` | 1,121 |

## Phenotypes / traits

| phenotype_id | distinct GSM | rows | numeric | categorical |
| --- | ---: | ---: | ---: | ---: |
| `tissue` | 169,869 | 232,417 | 0 | 232,417 |
| `sex` | 119,808 | 134,106 | 0 | 129,334 |
| `age` | 105,381 | 163,852 | 152,349 | 0 |
| `disease` | 97,732 | 104,082 | 0 | 97,152 |
| `bmi` | 33,652 | 45,846 | 6,358 | 0 |
| `cancer` | 15,324 | 15,324 | 0 | 13,447 |
| `blood` | 3,402 | 3,402 | 0 | 38 |
| `brain` | 1,997 | 1,997 | 0 | 1,997 |
| `ancestry` | 1,380 | 1,380 | 0 | 1,380 |

- Age (numeric): **96,031** GSM
- Sex: **119,808** GSM
- Tissue (mapped categorical in catalog): **169,869** GSM

### Top tissues

| Tissue | GSM |
| --- | ---: |
| whole blood | 46,001 |
| peripheral blood mononuclear cell | 8,425 |
| brain | 8,176 |
| leukocyte | 8,144 |
| breast | 5,366 |
| bone marrow | 3,960 |
| liver | 3,612 |
| cord blood | 3,464 |
| lung | 3,420 |
| peripheral blood | 3,398 |
| peripheral blood leukocyte | 3,244 |
| placenta | 2,554 |
| kidney | 2,519 |
| saliva | 2,426 |
| semen | 2,197 |
| meningioma tissue | 2,036 |
| colon | 1,928 |
| brain - glioblastoma tissue | 1,844 |
| CD14+ monocyte | 1,769 |
| brain - cerebellum | 1,610 |
| prostate gland | 1,575 |
| buccal epithelial cell | 1,556 |
| sperm | 1,528 |
| sarcoma | 1,504 |
| uterus - endometrium | 1,446 |

### Disease / cancer label status

**disease**

| label_status | GSM |
| --- | ---: |
| `case` | 80,133 |
| `control` | 10,960 |
| `unknown` | 6,930 |

**cancer**

| label_status | GSM |
| --- | ---: |
| `case` | 12,750 |
| `unknown` | 1,877 |
| `control` | 697 |

## Trait eligibility (training gates)

| phenotype | family | task | n | cases | controls | unknown | core | aux | external | exclusion |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- | --- | --- | --- |
| `age` | age | continuous | 8374 | None | None | 0 | yes | yes | no | nan |
| `bmi` | age | continuous | 879 | None | None | 7495 | no | yes | no | need ≥1000 samples, ≥5 studies, range across >1 study |
| `sex` | age | binary | 8271 | None | None | 103 | no | yes | no | sex is auxiliary biological / QC, not a core burden target |
| `tissue` | age | multiclass | 8374 | None | None | 0 | yes | yes | no | nan |
| `age` | ancestry | continuous | 1205 | None | None | 175 | yes | yes | no | nan |
| `ancestry` | ancestry | multiclass | 1380 | None | None | 0 | no | no | yes | ancestry is fairness / domain eval only |
| `bmi` | ancestry | continuous | 213 | None | None | 1167 | no | yes | no | need ≥1000 samples, ≥5 studies, range across >1 study |
| `sex` | ancestry | binary | 1380 | None | None | 0 | no | yes | no | sex is auxiliary biological / QC, not a core burden target |
| `tissue` | ancestry | multiclass | 1380 | None | None | 0 | no | yes | no | need ≥2 classes with ≥100 samples and ≥2 studies |
| `age` | blood | continuous | 1872 | None | None | 1530 | yes | yes | no | nan |
| `blood` | blood | multiclass | 38 | None | None | 3364 | no | no | no | fine-grained blood/brain outside single-study default core |
| `bmi` | blood | continuous | 44 | None | None | 3358 | no | no | no | need ≥1000 samples, ≥5 studies, range across >1 study |
| `sex` | blood | binary | 2478 | None | None | 924 | no | yes | no | sex is auxiliary biological / QC, not a core burden target |
| `tissue` | blood | multiclass | 3402 | None | None | 0 | yes | yes | no | nan |
| `age` | bmi | continuous | 2068 | None | None | 2 | yes | yes | no | nan |
| `bmi` | bmi | continuous | 2070 | None | None | 0 | yes | yes | no | nan |
| `sex` | bmi | binary | 2070 | None | None | 0 | no | yes | no | sex is auxiliary biological / QC, not a core burden target |
| `tissue` | bmi | multiclass | 2070 | None | None | 0 | yes | yes | no | nan |
| `age` | brain | continuous | 1738 | None | None | 259 | yes | yes | no | nan |
| `brain` | brain | multiclass | 1997 | None | None | 0 | no | yes | yes | fine-grained blood/brain outside single-study default core |
| `sex` | brain | binary | 1853 | None | None | 144 | no | yes | no | sex is auxiliary biological / QC, not a core burden target |
| `age` | cancer | continuous | 7717 | None | None | 2384 | yes | yes | no | nan |
| `bmi` | cancer | continuous | 1702 | None | None | 8399 | yes | yes | no | nan |
| `cancer` | cancer | binary_or_multilabel | 8224 | 8141 | 83 | 1877 | no | yes | yes | need ≥200 cases, ≥200 controls, ≥3 studies (unknown≠control) |
| `sex` | cancer | binary | 8780 | None | None | 1321 | no | yes | no | sex is auxiliary biological / QC, not a core burden target |
| `tissue` | cancer | multiclass | 10101 | None | None | 0 | yes | yes | no | nan |
| `age` | disease | continuous | 9373 | None | None | 2845 | yes | yes | no | nan |
| `bmi` | disease | continuous | 305 | None | None | 11913 | no | yes | no | need ≥1000 samples, ≥5 studies, range across >1 study |
| `disease` | disease | binary_or_multilabel | 5288 | 5288 | 0 | 6930 | no | yes | yes | need ≥200 cases, ≥200 controls, ≥3 studies (unknown≠control) |
| `sex` | disease | binary | 10968 | None | None | 1250 | no | yes | no | sex is auxiliary biological / QC, not a core burden target |
| `tissue` | disease | multiclass | 12218 | None | None | 0 | yes | yes | no | nan |
| `age` | ewas_datahub_repository | continuous | 66031 | None | None | 0 | yes | yes | no | nan |
| `disease` | ewas_datahub_repository | binary_or_multilabel | 73705 | 73705 | 0 | 0 | no | yes | yes | need ≥200 cases, ≥200 controls, ≥3 studies (unknown≠control) |
| `tissue` | ewas_datahub_repository | multiclass | 134089 | None | None | 0 | yes | yes | no | nan |
| `age` | geo_metadata_backfill | continuous | 49978 | None | None | 0 | yes | yes | no | nan |
| `cancer` | geo_metadata_backfill | binary_or_multilabel | 5223 | 4609 | 614 | 0 | yes | yes | yes | nan |
| `disease` | geo_metadata_backfill | binary_or_multilabel | 18159 | 7199 | 10960 | 0 | yes | yes | yes | nan |
| `sex` | geo_metadata_backfill | binary | 86263 | None | None | 0 | no | yes | no | sex is auxiliary biological / QC, not a core burden target |
| `tissue` | geo_metadata_backfill | multiclass | 52482 | None | None | 0 | yes | yes | no | nan |
| `age` | sex | continuous | 1617 | None | None | 1361 | yes | yes | no | nan |
| `bmi` | sex | continuous | 440 | None | None | 2538 | no | yes | no | need ≥1000 samples, ≥5 studies, range across >1 study |
| `sex` | sex | binary | 2978 | None | None | 0 | no | yes | no | sex is auxiliary biological / QC, not a core burden target |
| `tissue` | sex | multiclass | 2978 | None | None | 0 | yes | yes | no | nan |
| `age` | tissue | continuous | 2376 | None | None | 2947 | yes | yes | no | nan |
| `bmi` | tissue | continuous | 705 | None | None | 4618 | no | yes | no | need ≥1000 samples, ≥5 studies, range across >1 study |
| `sex` | tissue | binary | 4293 | None | None | 1030 | no | yes | no | sex is auxiliary biological / QC, not a core burden target |
| `tissue` | tissue | multiclass | 5323 | None | None | 0 | yes | yes | no | nan |

## Phenotype provenance (top source_family × phenotype)

| source_family | phenotype_id | GSM | rows |
| --- | --- | ---: | ---: |
| `ewas_datahub_repository` | `tissue` | 134,089 | 134,089 |
| `geo_metadata_backfill` | `sex` | 86,263 | 86,263 |
| `ewas_datahub_repository` | `disease` | 73,705 | 73,705 |
| `ewas_datahub_repository` | `age` | 66,031 | 66,031 |
| `geo_metadata_backfill` | `tissue` | 52,482 | 52,482 |
| `geo_metadata_backfill` | `age` | 49,978 | 49,978 |
| `geo_metadata_backfill` | `disease` | 18,159 | 18,159 |
| `disease` | `disease` | 12,218 | 12,218 |
| `disease` | `bmi` | 12,218 | 12,218 |
| `disease` | `sex` | 12,218 | 12,218 |
| `disease` | `tissue` | 12,218 | 12,218 |
| `disease` | `age` | 12,218 | 12,218 |
| `cancer` | `tissue` | 10,101 | 10,101 |
| `cancer` | `cancer` | 10,101 | 10,101 |
| `cancer` | `bmi` | 10,101 | 10,101 |
| `cancer` | `sex` | 10,101 | 10,101 |
| `cancer` | `age` | 10,101 | 10,101 |
| `age` | `tissue` | 8,374 | 8,374 |
| `age` | `bmi` | 8,374 | 8,374 |
| `age` | `sex` | 8,374 | 8,374 |
| `age` | `age` | 8,374 | 8,374 |
| `tissue` | `bmi` | 5,323 | 5,323 |
| `tissue` | `tissue` | 5,323 | 5,323 |
| `tissue` | `age` | 5,323 | 5,323 |
| `tissue` | `sex` | 5,323 | 5,323 |
| `geo_metadata_backfill` | `cancer` | 5,223 | 5,223 |
| `blood` | `tissue` | 3,402 | 3,402 |
| `blood` | `bmi` | 3,402 | 3,402 |
| `blood` | `sex` | 3,402 | 3,402 |
| `blood` | `age` | 3,402 | 3,402 |
| `blood` | `blood` | 3,402 | 3,402 |
| `sex` | `sex` | 2,978 | 2,978 |
| `sex` | `tissue` | 2,978 | 2,978 |
| `sex` | `bmi` | 2,978 | 2,978 |
| `sex` | `age` | 2,978 | 2,978 |
| `bmi` | `sex` | 2,070 | 2,070 |
| `bmi` | `tissue` | 2,070 | 2,070 |
| `bmi` | `age` | 2,070 | 2,070 |
| `bmi` | `bmi` | 2,070 | 2,070 |
| `brain` | `brain` | 1,997 | 1,997 |

## EWAS_db assay mirror (on-disk methylation profiles)

| Advertised | Local dirs | Dirs with GSM*.txt | GSM*.txt | Any *.txt |
| ---: | ---: | ---: | ---: | ---: |
| 1989 | 1989 | 1811 | 173,099 | 187,324 |

Empty-dir refill: `scripts/refill_ewas_db_empty_studies.sh` + `configs/data/ewas_db_refill_*.txt` (see `ewas_db_empty_refill_plan.json`).

## GEO metadata (labels / study context — not assay betas)

- Sample parquet: **170,338** GSM / **1707** primary studies; age **64,727**, sex **110,430**, tissue_map `{'unmapped': 100086, 'mapped': 70244, 'empty': 8}`
- Series: **1813** GSE; title **1718**; overall_design **271**

## Related reports

- [`census.md`](census.md)
- [`trait_eligibility.md`](trait_eligibility.md)
- [`geo_series_enrichment.md`](geo_series_enrichment.md)
- [`ewas_db_download_failures.md`](ewas_db_download_failures.md)

