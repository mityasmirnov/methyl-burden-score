# EWAS metadata structure

Generated at: `2026-08-06T14:48:43.718354+00:00`

Scope: small EWAS Atlas tables + unpacked DataHub `sample_*.txt` packs.
Large Atlas associations / probe annotations and matrix zips are out of scope.

## Parse recipes

| Source | Recipe |
|--------|--------|
| Atlas studies / cohorts | TSV, header row, tab-separated |
| Atlas trait×trait | TSV matrix; col0 = trait; remaining = traits |
| DataHub sample-info | R write.table (space + quotes); `read_r_style_table` |

## Family → primary phenotype column

| Family | Column |
|--------|--------|
| `age` | `age` |
| `ancestry` | `race` |
| `blood` | `cell_component` |
| `bmi` | `bmi` |
| `brain` | `tissue` |
| `cancer` | `disease` |
| `disease` | `disease` |
| `sex` | `sex` |
| `tissue` | `tissue` |

## Atlas small tables

### `studies`

- Path: `/data/projects/methyl-burden-score/data/raw/ewas_atlas/EWAS_Atlas_studies.tsv`
- Bytes: 153681
- Shape: 1902 × 5
- Parse: tab-separated TSV with header (latin-1); skip rare malformed rows
- Join keys: `study_ID`, `PMID`

| Column | Kind | Non-null | N unique | Notes |
|--------|------|---------:|---------:|-------|
| `study_ID` | string | 1.000 | 1902 |  |
| `trait` | string | 1.000 | 878 |  |
| `case_description` | string | 0.910 | 1291 |  |
| `control_description` | string | 0.503 | 453 |  |
| `PMID` | numeric | 1.000 | 1182 | min=20019873.0, max=39827095.0 |

Example rows:

```json
[
  {
    "study_ID": "ES00033",
    "trait": "body mass index (BMI)",
    "case_description": null,
    "control_description": null,
    "PMID": "24630777"
  },
  {
    "study_ID": "ES00034",
    "trait": "body mass index (BMI)",
    "case_description": null,
    "control_description": null,
    "PMID": "24630777"
  }
]
```

### `cohorts`

- Path: `/data/projects/methyl-burden-score/data/raw/ewas_atlas/EWAS_Atlas_cohorts.tsv`
- Bytes: 397573
- Shape: 3982 × 15
- Parse: tab-separated TSV with header (latin-1); skip rare malformed rows
- Malformed rows skipped: 1
- Join keys: `study_ID`, `cohort_ID`

| Column | Kind | Non-null | N unique | Notes |
|--------|------|---------:|---------:|-------|
| `cohort_ID` | numeric | 1.000 | 3982 | min=52.0, max=4563.0 |
| `study_ID` | string | 1.000 | 1788 |  |
| `stage` | categorical | 1.000 | 2 | Initial (3543), Replication (439) |
| `platform` | categorical | 0.997 | 3 | 450K (3134), 850K (710), 27K (126) |
| `sample_size` | numeric | 1.000 | 891 | min=1.0, max=22774.0 |
| `male_percentage` | string | 0.774 | 257 |  |
| `min_age` | categorical | 0.325 | 151 | 0 (432), 14 (56), 18 (34) |
| `max_age` | categorical | 0.325 | 173 | 0 (420), 94 (52), 81 (33) |
| `mean_age` | string | 0.753 | 736 |  |
| `sd_age` | string | 0.609 | 347 |  |
| `tissue` | string | 0.999 | 250 |  |
| `cohort_name` | string | 0.490 | 487 |  |
| `full_name` | string | 0.285 | 288 |  |
| `description` | string | 0.455 | 810 |  |
| `ancestry` | categorical | 0.997 | 101 | Not reported (1719), European (1262), East Asian (189) |

Example rows:

```json
[
  {
    "cohort_ID": "52",
    "study_ID": "ES00033",
    "stage": "Initial",
    "platform": "450K",
    "sample_size": "239",
    "male_percentage": "0.85",
    "min_age": null,
    "max_age": null,
    "mean_age": "55.2",
    "sd_age": "6.8",
    "tissue": "whole blood",
    "cohort_name": null,
    "full_name": "Cardiogenics Consortium",
    "description": "myocardial infarction",
    "ancestry": "European"
  },
  {
    "cohort_ID": "53",
    "study_ID": "ES00033",
    "stage": "Initial",
    "platform": "450K",
    "sample_size": "220",
    "male_percentage": "0.57",
    "min_age": null,
    "max_age": null,
    "mean_age": "55.2",
    "sd_age": "6.8",
    "tissue": "whole blood",
    "cohort_name": null,
    "full_name": "Cardiogenics Consortium",
    "description": "Healthy blood donors",
    "ancestry": "European"
  }
]
```

### `trait_trait_logP`

- Path: `/data/projects/methyl-burden-score/data/raw/ewas_atlas/EWAS_trait_trait_logP.txt`
- Bytes: 986305
- Shape: 372 × 373
- Parse: tab-separated latin-1; first column = trait name; remaining = square logP matrix
- Join keys: `trait (row/col labels; not Atlas study_ID)`
- Traits: 372; square=True
- Value sample range: [-0.0, 400.0]


## DataHub sample-info packs

### Family `age`

- Path: `/data/projects/methyl-burden-score/reports/inspection/ewas_datahub_samples/sample_age_methylation_v1/sample_age.txt`
- Bytes: 2287958
- Shape: 8374 × 67
- Join keys: `sample_id`, `project_id`
- Primary phenotype: `age`
- Primary non-null=1.0, n_unique=1734, kind=numeric
- Numeric range: [0.0, 114.0] mean=49.9758344284213
- N sample_id / project_id: 8374 / 143

### Family `ancestry`

- Path: `/data/projects/methyl-burden-score/reports/inspection/ewas_datahub_samples/sample_ancestry_category_methylation_v1/sample_race.txt`
- Bytes: 224431
- Shape: 1380 × 23
- Join keys: `sample_id`, `project_id`
- Primary phenotype: `race`
- Primary non-null=1.0, n_unique=16, kind=categorical
- Top values: African American (237), Hispanic - Mexican (227), Hispanic (168), Chinese (167), white (159)
- N sample_id / project_id: 1380 / 21

### Family `blood`

- Path: `/data/projects/methyl-burden-score/reports/inspection/ewas_datahub_samples/sample_blood_methylation_v1/sample_blood.txt`
- Bytes: 933472
- Shape: 3402 × 69
- Join keys: `sample_id`, `project_id`
- Primary phenotype: `cell_component`
- Primary non-null=0.0112, n_unique=38, kind=categorical
- Top values: CD8 T cell:0, CD4 T cell:0.0583, NK cell:0.0829, B cell:0.01947, monocyte:0.08806, granulocyte:0.76116 (1), CD8 T cell:0.18733, CD4 T cell:0.05062, NK cell:0.03571, B cell:0.04623, monocyte:0.1001, granulocyte:0.59217 (1), CD8 T cell:0.04315, CD4 T cell:0.09799, NK cell:0.03696, B cell:0.02177, monocyte:0.09522, granulocyte:0.70113 (1), CD8 T cell:0.13038, CD4 T cell:0.24836, NK cell:0.13716, B cell:0.07829, monocyte:0.09094, granulocyte:0.34153 (1), CD8 T cell:0.06984, CD4 T cell:0.19197, NK cell:0.06604, B cell:0.07694, monocyte:0.10896, granulocyte:0.51312 (1)
- N sample_id / project_id: 3402 / 161

### Family `bmi`

- Path: `/data/projects/methyl-burden-score/reports/inspection/ewas_datahub_samples/sample_bmi_methylation_v1/sample_bmi.txt`
- Bytes: 483837
- Shape: 2070 × 38
- Join keys: `sample_id`, `project_id`
- Primary phenotype: `bmi`
- Primary non-null=1.0, n_unique=1451, kind=numeric
- Numeric range: [6.07, 73.617] mean=30.873200966183575
- N sample_id / project_id: 2070 / 25

### Family `brain`

- Path: `/data/projects/methyl-burden-score/reports/inspection/ewas_datahub_samples/sample_brain_methylation_v1/sample_brain.txt`
- Bytes: 257892
- Shape: 1997 × 19
- Join keys: `sample_id`, `project_id`
- Primary phenotype: `tissue`
- Primary non-null=1.0, n_unique=41, kind=categorical
- Top values: brain - cerebellum (300), brain - dorsolateral prefrontal cortex (245), brain - superior temporal gyrus (179), brain - frontal lobe (155), brain - frontal cortex (133)
- N sample_id / project_id: 1997 / 40

### Family `cancer`

- Path: `/data/projects/methyl-burden-score/reports/inspection/ewas_datahub_samples/sample_cancer_methylation_v1/sample_cancer.txt`
- Bytes: 4005665
- Shape: 10841 × 73
- Join keys: `sample_id`, `project_id`
- Primary phenotype: `disease`
- Primary non-null=0.7999, n_unique=43, kind=categorical
- Top values: acute myeloid leukemia (657), prostate cancer (613), breast cancer (533), glioma (481), head and neck squamous-cell carcinoma (470)
- N sample_id / project_id: 10101 / 225

### Family `disease`

- Path: `/data/projects/methyl-burden-score/reports/inspection/ewas_datahub_samples/sample_disease_methylation_v1/sample_disease.txt`
- Bytes: 5491933
- Shape: 14501 × 100
- Join keys: `sample_id`, `project_id`
- Primary phenotype: `disease`
- Primary non-null=0.3664, n_unique=30, kind=categorical
- Top values: Alzheimer's disease (945), schizophrenia (536), systemic lupus erythematosus (341), Parkinson's disease (333), rheumatoid arthritis (250)
- N sample_id / project_id: 12218 / 209

### Family `sex`

- Path: `/data/projects/methyl-burden-score/reports/inspection/ewas_datahub_samples/sample_sex_methylation_v1/sample_sex.txt`
- Bytes: 991227
- Shape: 2978 × 84
- Join keys: `sample_id`, `project_id`
- Primary phenotype: `sex`
- Primary non-null=1.0, n_unique=2, kind=categorical
- Top values: F (1543), M (1435)
- N sample_id / project_id: 2978 / 161

### Family `tissue`

- Path: `/data/projects/methyl-burden-score/reports/inspection/ewas_datahub_samples/sample_tissue_methylation_v1/sample_tissue.txt`
- Bytes: 2069568
- Shape: 5323 × 104
- Join keys: `sample_id`, `project_id`
- Primary phenotype: `tissue`
- Primary non-null=1.0, n_unique=72, kind=categorical
- Top values: whole blood (300), cord blood (300), saliva (300), liver (300), breast (300)
- N sample_id / project_id: 5323 / 258

## Cross-pack column / ID overlap

- Packs profiled: 9
- Shared columns (10): `age`, `genotype`, `infection`, `platform`, `project_id`, `race`, `sample_id`, `sample_type`, `sex`, `tissue`

Pairwise `sample_id` overlap:

| A | B | Shared | N_A | N_B |
|---|---|-------:|----:|----:|
| `age` | `ancestry` | 220 | 8374 | 1380 |
| `age` | `blood` | 955 | 8374 | 3402 |
| `age` | `bmi` | 854 | 8374 | 2070 |
| `age` | `brain` | 1071 | 8374 | 1997 |
| `age` | `cancer` | 897 | 8374 | 10101 |
| `age` | `disease` | 2879 | 8374 | 12218 |
| `age` | `sex` | 570 | 8374 | 2978 |
| `age` | `tissue` | 931 | 8374 | 5323 |
| `ancestry` | `blood` | 24 | 1380 | 3402 |
| `ancestry` | `bmi` | 213 | 1380 | 2070 |
| `ancestry` | `brain` | 0 | 1380 | 1997 |
| `ancestry` | `cancer` | 42 | 1380 | 10101 |
| `ancestry` | `disease` | 229 | 1380 | 12218 |
| `ancestry` | `sex` | 14 | 1380 | 2978 |
| `ancestry` | `tissue` | 21 | 1380 | 5323 |
| `blood` | `bmi` | 44 | 3402 | 2070 |
| `blood` | `brain` | 0 | 3402 | 1997 |
| `blood` | `cancer` | 83 | 3402 | 10101 |
| `blood` | `disease` | 1033 | 3402 | 12218 |
| `blood` | `sex` | 70 | 3402 | 2978 |
| `blood` | `tissue` | 91 | 3402 | 5323 |
| `bmi` | `brain` | 0 | 2070 | 1997 |
| `bmi` | `cancer` | 227 | 2070 | 10101 |
| `bmi` | `disease` | 233 | 2070 | 12218 |
| `bmi` | `sex` | 372 | 2070 | 2978 |
| `bmi` | `tissue` | 612 | 2070 | 5323 |
| `brain` | `cancer` | 84 | 1997 | 10101 |
| `brain` | `disease` | 1041 | 1997 | 12218 |
| `brain` | `sex` | 0 | 1997 | 2978 |
| `brain` | `tissue` | 0 | 1997 | 5323 |
| `cancer` | `disease` | 287 | 10101 | 12218 |
| `cancer` | `sex` | 928 | 10101 | 2978 |
| `cancer` | `tissue` | 1422 | 10101 | 5323 |
| `disease` | `sex` | 411 | 12218 | 2978 |
| `disease` | `tissue` | 523 | 12218 | 5323 |
| `sex` | `tissue` | 1965 | 2978 | 5323 |

Pairwise `project_id` overlap:

| A | B | Shared | N_A | N_B |
|---|---|-------:|----:|----:|
| `age` | `ancestry` | 13 | 143 | 21 |
| `age` | `blood` | 62 | 143 | 161 |
| `age` | `bmi` | 16 | 143 | 25 |
| `age` | `brain` | 23 | 143 | 40 |
| `age` | `cancer` | 87 | 143 | 225 |
| `age` | `disease` | 111 | 143 | 209 |
| `age` | `sex` | 59 | 143 | 161 |
| `age` | `tissue` | 80 | 143 | 258 |
| `ancestry` | `blood` | 12 | 21 | 161 |
| `ancestry` | `bmi` | 2 | 21 | 25 |
| `ancestry` | `brain` | 0 | 21 | 40 |
| `ancestry` | `cancer` | 13 | 21 | 225 |
| `ancestry` | `disease` | 16 | 21 | 209 |
| `ancestry` | `sex` | 12 | 21 | 161 |
| `ancestry` | `tissue` | 14 | 21 | 258 |
| `blood` | `bmi` | 4 | 161 | 25 |
| `blood` | `brain` | 4 | 161 | 40 |
| `blood` | `cancer` | 70 | 161 | 225 |
| `blood` | `disease` | 108 | 161 | 209 |
| `blood` | `sex` | 52 | 161 | 161 |
| `blood` | `tissue` | 72 | 161 | 258 |
| `bmi` | `brain` | 0 | 25 | 40 |
| `bmi` | `cancer` | 11 | 25 | 225 |
| `bmi` | `disease` | 7 | 25 | 209 |
| `bmi` | `sex` | 15 | 25 | 161 |
| `bmi` | `tissue` | 22 | 25 | 258 |
| `brain` | `cancer` | 16 | 40 | 225 |
| `brain` | `disease` | 29 | 40 | 209 |
| `brain` | `sex` | 4 | 40 | 161 |
| `brain` | `tissue` | 10 | 40 | 258 |
| `cancer` | `disease` | 86 | 225 | 209 |
| `cancer` | `sex` | 89 | 225 | 161 |
| `cancer` | `tissue` | 122 | 225 | 258 |
| `disease` | `sex` | 82 | 209 | 161 |
| `disease` | `tissue` | 108 | 209 | 258 |
| `sex` | `tissue` | 153 | 161 | 258 |

## Atlas study_ID vs Hub project_id

- Atlas study_ID count: 1902
- Hub project_id count (union): 470
- Hub GSE-like: 435
- Exact string equals: 0
- Note: Atlas study_ID (ES*) and Hub project_id (usually GSE*) are different namespaces; do not join on raw equality. Use PMID / curated maps when needed.

## Related

- Durable contracts: `docs/EWAS_METADATA.md`
- Export path: `mbs.registry.sample_info`
