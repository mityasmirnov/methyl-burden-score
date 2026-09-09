# Data population inventory

- Generated: `2026-09-09T09:43:40Z`
- Release: `deepmat-data-v1`

## Catalog scale

| Samples | Studies | Phenotype rows | Assay files |
| ---: | ---: | ---: | ---: |
| **173,076** | **1,763** | **673,786** | **170,641** |

## Assays / platforms

Platform is study-level (`study.platform_id`); counts are samples under those studies.

| Platform | Samples |
| --- | ---: |
| `HM450` | 105,639 |
| `EPIC` | 55,901 |
| `(null)` | 10,415 |
| `EPICv2` | 1,121 |

## Phenotypes / traits

| phenotype_id | distinct GSM | rows | numeric | categorical |
| --- | ---: | ---: | ---: | ---: |
| `tissue` | 157,909 | 218,900 | 0 | 218,900 |
| `sex` | 117,671 | 131,969 | 0 | 127,197 |
| `age` | 100,260 | 157,911 | 146,408 | 0 |
| `disease` | 90,877 | 97,134 | 0 | 90,204 |
| `bmi` | 33,652 | 45,846 | 6,358 | 0 |
| `cancer` | 15,247 | 15,247 | 0 | 13,370 |
| `blood` | 3,402 | 3,402 | 0 | 38 |
| `brain` | 1,997 | 1,997 | 0 | 1,997 |
| `ancestry` | 1,380 | 1,380 | 0 | 1,380 |

- Age (numeric): **90,910** GSM
- Sex: **117,671** GSM
- Tissue (mapped categorical in catalog): **157,909** GSM

### Top tissues

| Tissue | GSM |
| --- | ---: |
| whole blood | 41,934 |
| peripheral blood mononuclear cell | 8,345 |
| leukocyte | 8,117 |
| brain | 8,045 |
| breast | 5,362 |
| peripheral blood | 3,308 |
| lung | 3,274 |
| cord blood | 3,205 |
| liver | 3,105 |
| peripheral blood leukocyte | 3,053 |
| placenta | 2,541 |
| kidney | 2,485 |
| saliva | 2,365 |
| semen | 2,194 |
| bone marrow | 2,075 |
| meningioma tissue | 2,036 |
| colon | 1,880 |
| brain - glioblastoma tissue | 1,765 |
| CD14+ monocyte | 1,760 |
| buccal epithelial cell | 1,556 |
| sperm | 1,528 |
| sarcoma | 1,504 |
| prostate gland | 1,464 |
| uterus - endometrium | 1,420 |
| buccal | 1,342 |

### Disease / cancer label status

**disease**

| label_status | GSM |
| --- | ---: |
| `case` | 73,711 |
| `control` | 10,441 |
| `unknown` | 6,930 |

**cancer**

| label_status | GSM |
| --- | ---: |
| `case` | 12,682 |
| `unknown` | 1,877 |
| `control` | 688 |

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
| `age` | ewas_datahub_repository | continuous | 61078 | None | None | 0 | yes | yes | no | nan |
| `disease` | ewas_datahub_repository | binary_or_multilabel | 67283 | 67283 | 0 | 0 | no | yes | yes | need ≥200 cases, ≥200 controls, ≥3 studies (unknown≠control) |
| `tissue` | ewas_datahub_repository | multiclass | 122130 | None | None | 0 | yes | yes | no | nan |
| `age` | geo_metadata_backfill | continuous | 48990 | None | None | 0 | yes | yes | no | nan |
| `cancer` | geo_metadata_backfill | binary_or_multilabel | 5146 | 4541 | 605 | 0 | yes | yes | yes | nan |
| `disease` | geo_metadata_backfill | binary_or_multilabel | 17633 | 7192 | 10441 | 0 | yes | yes | yes | nan |
| `sex` | geo_metadata_backfill | binary | 84126 | None | None | 0 | no | yes | no | sex is auxiliary biological / QC, not a core burden target |
| `tissue` | geo_metadata_backfill | multiclass | 50924 | None | None | 0 | yes | yes | no | nan |
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
| `ewas_datahub_repository` | `tissue` | 122,130 | 122,130 |
| `geo_metadata_backfill` | `sex` | 84,126 | 84,126 |
| `ewas_datahub_repository` | `disease` | 67,283 | 67,283 |
| `ewas_datahub_repository` | `age` | 61,078 | 61,078 |
| `geo_metadata_backfill` | `tissue` | 50,924 | 50,924 |
| `geo_metadata_backfill` | `age` | 48,990 | 48,990 |
| `geo_metadata_backfill` | `disease` | 17,633 | 17,633 |
| `disease` | `disease` | 12,218 | 12,218 |
| `disease` | `sex` | 12,218 | 12,218 |
| `disease` | `age` | 12,218 | 12,218 |
| `disease` | `bmi` | 12,218 | 12,218 |
| `disease` | `tissue` | 12,218 | 12,218 |
| `cancer` | `cancer` | 10,101 | 10,101 |
| `cancer` | `age` | 10,101 | 10,101 |
| `cancer` | `tissue` | 10,101 | 10,101 |
| `cancer` | `sex` | 10,101 | 10,101 |
| `cancer` | `bmi` | 10,101 | 10,101 |
| `age` | `bmi` | 8,374 | 8,374 |
| `age` | `tissue` | 8,374 | 8,374 |
| `age` | `age` | 8,374 | 8,374 |
| `age` | `sex` | 8,374 | 8,374 |
| `tissue` | `sex` | 5,323 | 5,323 |
| `tissue` | `bmi` | 5,323 | 5,323 |
| `tissue` | `tissue` | 5,323 | 5,323 |
| `tissue` | `age` | 5,323 | 5,323 |
| `geo_metadata_backfill` | `cancer` | 5,146 | 5,146 |
| `blood` | `age` | 3,402 | 3,402 |
| `blood` | `tissue` | 3,402 | 3,402 |
| `blood` | `bmi` | 3,402 | 3,402 |
| `blood` | `sex` | 3,402 | 3,402 |
| `blood` | `blood` | 3,402 | 3,402 |
| `sex` | `tissue` | 2,978 | 2,978 |
| `sex` | `bmi` | 2,978 | 2,978 |
| `sex` | `age` | 2,978 | 2,978 |
| `sex` | `sex` | 2,978 | 2,978 |
| `bmi` | `age` | 2,070 | 2,070 |
| `bmi` | `bmi` | 2,070 | 2,070 |
| `bmi` | `tissue` | 2,070 | 2,070 |
| `bmi` | `sex` | 2,070 | 2,070 |
| `brain` | `sex` | 1,997 | 1,997 |

## EWAS_db assay mirror (on-disk methylation profiles)

| Advertised | Local dirs | Dirs with GSM*.txt | GSM*.txt | Any *.txt |
| ---: | ---: | ---: | ---: | ---: |
| 1989 | 1989 | 1655 | 157,279 | 170,747 |

Gap fill (running): empty dirs with remote content — **162** GSE (~15.9k GSM files) + **49** non-GSM namespaces (~14.2k TCGA/ArrayExpress/TARGET/ENCODE/… files). See `ewas_db_empty_refill_plan.json`, `scripts/refill_ewas_db_empty_studies.sh`.

## GEO metadata (labels / study context — not assay betas)

- Sample parquet: **170,338** GSM / **1707** primary studies; age **64,727**, sex **110,430**, tissue_map `{'unmapped': 100086, 'mapped': 70244, 'empty': 8}`
- Series: **1718** GSE with title **1718**, overall_design **271**

## Related reports

- [`census.md`](census.md)
- [`trait_eligibility.md`](trait_eligibility.md)
- [`geo_series_enrichment.md`](geo_series_enrichment.md)
- [`ewas_db_download_failures.md`](ewas_db_download_failures.md)
- [`geo_backfill_remain/analysis.md`](geo_backfill_remain/analysis.md)

