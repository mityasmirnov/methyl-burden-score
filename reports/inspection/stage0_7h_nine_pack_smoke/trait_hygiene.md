# Milestone 10c — trait expansion hygiene

Policy + census only. **Does not** rewrite phenotype parquet or launch trait GPU.
Nine-pack samples: **34234**.

## Locked policy

1. **Disease/cancer case≠control:** use Hub `sample_type` mapped through `mbs.datahub_census.SAMPLE_TYPE_CASE_CONTROL` (`control`→control, `disease tissue`→case, `adjacent normal`→adjacent_normal). Train masked heads only on rows with `label_status ∈ {case, control}`. Pack-membership `disease_mask`/`cancer_mask` false means **unknown / not in pack**, **not** control. Diagnosis strings in `phenotype_value` remain multilabel targets.
2. **Blood/brain:** **defer** trait heads. Packs exist but nine-pack `blood_mask`/`brain_mask` are hardcoded false; sidecars label every row `sample_type=control` (pack catalogue, not case/control).
3. **Tissue:** do **not** expand the tissue CE head — only `whole blood` ≥1k among 64 labels. Any collapse map is a later ADR/config.

## Disease / cancer census (sample_type → label_status)

| Family | pack samples | case | control | adjacent_normal | ≥1k case+control? | nine-pack mask_true |
|--------|-------------:|-----:|--------:|----------------:|:-----------------:|--------------------:|
| disease | 12218 | 5264 | 6930 | 24 | yes | 12218 |
| cancer | 10101 | 7157 | 1920 | 1024 | yes | 10101 |

### sample_type raw counts

- disease: `{'control': 6930, 'disease tissue': 5264, 'adjacent normal': 24}`
- cancer: `{'disease tissue': 7157, 'control': 1920, 'adjacent normal': 1024}`

## Blood / brain

- blood pack n=3402; nine-pack `blood_mask` true=0; sample_type={'control': 3402} → **defer_head**
- brain pack n=1997; nine-pack `brain_mask` true=0; sample_type={'control': 1997} → **defer_head**

## Tissue collapse

- labeled=13457; labels=64
- ≥1k: `{'whole blood': 2070}`
- top5: `{'whole blood': 2070, 'leukocyte': 934, 'CD14+ monocyte': 724, 'saliva': 629, 'breast': 595}`
- decision: **do_not_expand_tissue_ce**

## Next (out of this slice)

- Implement policy in a new phenotype table / longform label_status columns (not now).
- After 10a+10b review, only then consider disease/cancer GPU with explicit case/control masks.
