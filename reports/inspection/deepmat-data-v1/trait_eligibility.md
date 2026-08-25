# Trait eligibility (deepmat-data-v1)

- Generated: `2026-08-25T11:21:35Z`

| Family | Phenotype | Task | n | studies | core | aux | ext | Reason |
| --- | --- | --- | ---: | ---: | --- | --- | --- | --- |
| `age` | `age` | continuous | 2 | 1 | False | False | False | need ≥1000 samples, ≥5 studies, range across >1 study |
| `age` | `sex` | binary | 2 | 1 | False | False | False | sex is auxiliary biological / QC, not a core burden target |
| `disease` | `disease` | binary_or_multilabel | 2 | 2 | False | False | False | need ≥200 cases, ≥200 controls, ≥3 studies (unknown≠control) |
| `tissue` | `sex` | binary | 2 | 2 | False | False | False | sex is auxiliary biological / QC, not a core burden target |
| `tissue` | `tissue` | multiclass | 2 | 2 | False | False | False | need ≥2 classes with ≥100 samples and ≥2 studies |

Disease/cancer: unknown labels are never treated as controls (ADR / DATA_CONTRACT).

