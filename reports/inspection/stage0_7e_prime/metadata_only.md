# Milestone 7E′ metadata-only ceiling

- Generated: `2026-08-25T12:35:49Z`
- Protocol: `fit_train_score_holdout`

Metadata-only confounding ceiling (study/platform/tissue one-hots). Sidecar for 7E′; 7E bake-off report must include the ATS ceiling. Not a replacement for neural training.

## ats_freeze_hub-age-tissue-sex-full-auto-v1

- split_id: `hub-age-tissue-sex-full-auto-v1`
- n_train / val / test: **9489** / **2074** / **1985**

### validation

- **age**: mae=20.02326631901986, rmse=23.612006842624922, r2=-0.14348288129355669, pearson_r=-0.12475298419394946, spearman_r=-0.1440155050963179
- **tissue**: macro_f1=0.3805463136410658, balanced_accuracy=0.3829787234042553
- **sex**: macro_f1=0.6169340237543454, balanced_accuracy=0.6243231820492511

### external_test

- **age**: mae=20.938508267267082, rmse=24.668475179571658, r2=0.08339392660345513, pearson_r=0.33419015028535026, spearman_r=0.09880672470618648
- **tissue**: macro_f1=0.3132714235335606, balanced_accuracy=0.3191489361702128
- **sex**: macro_f1=0.5317343186265872, balanced_accuracy=0.5479779739560879

## hub_nine_pack_full-auto-v1

- split_id: `hub-nine-pack-full-auto-v1`
- n_train / val / test: **24010** / **5576** / **4648**

### validation

- **age**: mae=19.25737525814981, rmse=22.55703532964939, r2=0.052743933684797684, pearson_r=0.251199866495667, spearman_r=0.1520700287601159
- **tissue**: macro_f1=0.47492458521870284, balanced_accuracy=0.484375
- **sex**: macro_f1=0.6298759606443736, balanced_accuracy=0.6311367676662047

### external_test

- **age**: mae=22.33882973737262, rmse=25.694626445014716, r2=0.052221651026038374, pearson_r=0.549792372138225, spearman_r=0.5627690590696226
- **tissue**: macro_f1=0.39748251104565535, balanced_accuracy=0.40625
- **sex**: macro_f1=0.4503280835106299, balanced_accuracy=0.451638457256758

## Virtual Hub store

- `matrix-hub-nine-pack-virtual-v1` (route + indices; no dense Zarr)
- Phenotype table: `sample_phenotype_table_hub_nine_pack_v1.parquet`
- Blood `cell_component` is **not** a pack-wide head

