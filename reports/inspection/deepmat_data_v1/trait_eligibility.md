# Trait eligibility (deepmat-data-v1)

- Generated: `2026-09-03T09:22:29Z`

| Family | Phenotype | Task | n | studies | core | aux | ext | Reason |
| --- | --- | --- | ---: | ---: | --- | --- | --- | --- |
| `age` | `age` | continuous | 8374 | 143 | True | True | False |  |
| `age` | `bmi` | continuous | 879 | 17 | False | True | False | need ≥1000 samples, ≥5 studies, range across >1 study |
| `age` | `sex` | binary | 8271 | 140 | False | True | False | sex is auxiliary biological / QC, not a core burden target |
| `age` | `tissue` | multiclass | 8374 | 143 | True | True | False |  |
| `ancestry` | `age` | continuous | 1205 | 17 | True | True | False |  |
| `ancestry` | `ancestry` | multiclass | 1380 | 21 | False | False | True | ancestry is fairness / domain eval only |
| `ancestry` | `bmi` | continuous | 213 | 2 | False | True | False | need ≥1000 samples, ≥5 studies, range across >1 study |
| `ancestry` | `sex` | binary | 1380 | 21 | False | True | False | sex is auxiliary biological / QC, not a core burden target |
| `ancestry` | `tissue` | multiclass | 1380 | 21 | False | True | False | need ≥2 classes with ≥100 samples and ≥2 studies |
| `blood` | `age` | continuous | 1872 | 74 | True | True | False |  |
| `blood` | `blood` | multiclass | 38 | 2 | False | False | False | fine-grained blood/brain outside single-study default core |
| `blood` | `bmi` | continuous | 44 | 4 | False | False | False | need ≥1000 samples, ≥5 studies, range across >1 study |
| `blood` | `sex` | binary | 2478 | 124 | False | True | False | sex is auxiliary biological / QC, not a core burden target |
| `blood` | `tissue` | multiclass | 3402 | 161 | True | True | False |  |
| `bmi` | `age` | continuous | 2068 | 25 | True | True | False |  |
| `bmi` | `bmi` | continuous | 2070 | 25 | True | True | False |  |
| `bmi` | `sex` | binary | 2070 | 25 | False | True | False | sex is auxiliary biological / QC, not a core burden target |
| `bmi` | `tissue` | multiclass | 2070 | 25 | True | True | False |  |
| `brain` | `age` | continuous | 1738 | 30 | True | True | False |  |
| `brain` | `brain` | multiclass | 1997 | 40 | False | True | True | fine-grained blood/brain outside single-study default core |
| `brain` | `sex` | binary | 1853 | 35 | False | True | False | sex is auxiliary biological / QC, not a core burden target |
| `cancer` | `age` | continuous | 7717 | 137 | True | True | False |  |
| `cancer` | `bmi` | continuous | 1702 | 22 | True | True | False |  |
| `cancer` | `cancer` | binary_or_multilabel | 8224 | 101 | False | True | True | need ≥200 cases, ≥200 controls, ≥3 studies (unknown≠control) |
| `cancer` | `sex` | binary | 8780 | 180 | False | True | False | sex is auxiliary biological / QC, not a core burden target |
| `cancer` | `tissue` | multiclass | 10101 | 225 | True | True | False |  |
| `disease` | `age` | continuous | 9373 | 133 | True | True | False |  |
| `disease` | `bmi` | continuous | 305 | 7 | False | True | False | need ≥1000 samples, ≥5 studies, range across >1 study |
| `disease` | `disease` | binary_or_multilabel | 5288 | 76 | False | True | True | need ≥200 cases, ≥200 controls, ≥3 studies (unknown≠control) |
| `disease` | `sex` | binary | 10968 | 188 | False | True | False | sex is auxiliary biological / QC, not a core burden target |
| `disease` | `tissue` | multiclass | 12218 | 209 | True | True | False |  |
| `geo_metadata_backfill` | `age` | continuous | 7652 | 5 | True | True | False |  |
| `geo_metadata_backfill` | `disease` | binary_or_multilabel | 2058 | 2 | False | False | False | need ≥200 cases, ≥200 controls, ≥3 studies (unknown≠control) |
| `geo_metadata_backfill` | `sex` | binary | 12729 | 8 | False | True | False | sex is auxiliary biological / QC, not a core burden target |
| `geo_metadata_backfill` | `tissue` | multiclass | 10528 | 7 | True | True | False |  |
| `sex` | `age` | continuous | 1617 | 95 | True | True | False |  |
| `sex` | `bmi` | continuous | 440 | 21 | False | True | False | need ≥1000 samples, ≥5 studies, range across >1 study |
| `sex` | `sex` | binary | 2978 | 161 | False | True | False | sex is auxiliary biological / QC, not a core burden target |
| `sex` | `tissue` | multiclass | 2978 | 161 | True | True | False |  |
| `tissue` | `age` | continuous | 2376 | 126 | True | True | False |  |
| `tissue` | `bmi` | continuous | 705 | 28 | False | True | False | need ≥1000 samples, ≥5 studies, range across >1 study |
| `tissue` | `sex` | binary | 4293 | 207 | False | True | False | sex is auxiliary biological / QC, not a core burden target |
| `tissue` | `tissue` | multiclass | 5323 | 258 | True | True | False |  |

Disease/cancer: unknown labels are never treated as controls (ADR / DATA_CONTRACT).

