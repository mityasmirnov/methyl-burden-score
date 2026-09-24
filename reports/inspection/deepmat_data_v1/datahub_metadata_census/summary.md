# DataHub metadata census

## Official vs local census Parquet

| Metric | Official | Local |
| --- | ---: | ---: |
| Samples | 180317 | 180317 |
| Tissues/Cells | 1413 | 1415 |
| Diseases | 842 | 847 |
| Fields | 296 | 269 |

## Platform split (local census)

- `450K`: 105085
- `850K`: 74187
- `935K`: 1045

## Catalog coverage

- Overlap census ∩ catalog: **169367**
- In census, not yet in catalog (missing betas): **10950**
- In catalog, not in census: **16484**
- Hub baseline samples present in census: **34234**

## Merge into refresh

- Enabled stats: `{"n_census_rows": 180317, "n_samples_matched": 169367, "n_platform_filled": 169367, "n_tissue_filled": 134671, "n_sample_type_filled": 155290, "n_phenotype_rows_added": 273825, "n_studies_platform_set": 1002, "enabled": true, "parquet_path": "/data/projects/methyl-burden-score/data/canonical/phenotypes/ewas_datahub_sample_census.parquet"}`

Hub packs remain training phenotype SoT; DataHub census fills platform/tissue/sample_type and `metadata_json.datahub` for stratification.
