# Residual catalog / DataHub / platform gaps

- Generated: `2026-09-07T14:33:00Z` (post-retry refresh)

## EWAS_db retry

- Clean retry finished; failures audit **`still_missing=0`** (historical WARN lines remain in the download log).
- On-disk EWAS_db: **1 695** studies with `.txt`, **170 641** sample files (**157 240** GSM-named).
- Release EWAS_db: `{"advertised_n": 1989, "mirror_complete": false, "n_local_gsm": 170641, "n_local_studies": 1695, "refresh_is_idempotent": true, "remote_index_fetched": false, "total_bytes": 1793230765292}`

## Catalog vs DataHub census

| Metric | N |
| --- | ---: |
| Official DataHub samples | 180317 |
| Local census Parquet | 180317 |
| Catalog samples | 173076 |
| Overlap census ∩ catalog | 157402 |
| In census, not in catalog (missing betas) | 22915 |
| In catalog, not in census | 15674 |

## Study `platform_id`

- Filled: **1456 / 1763** (82.6%)
- `HM450`: 908
- `EPIC`: 530
- `nan`: 307
- `EPICv2`: 18

- Still null: **307** studies (13321 samples)
  - by kind: `{'GSE': 304, 'ArrayExpress': 2, 'TCGA/CPTAC': 1}`
  - reasons: `{'mixed_census_or_geo': 43, 'single_unapplied_zero_samples': 3, 'not_in_census_or_geo': 261}`

## Hub lane flags

- `{'hub_only': 3900, 'ewas_db_only': 138842, 'both': 30334}`

- Registry `sample_count: null` remaining: **0**

## Notes

- Hub pack phenotypes remain training SoT.
- Mixed-platform studies correctly keep study.platform_id null.
- single_unapplied_zero_samples = inventory studies with no catalog samples but a single census/GEO platform.
- GEO priority gap list filled platforms where SOFT had a single methylation GPL.
