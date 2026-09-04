# Seed-panel audit (fold-0)

- panel_hash: `ef6cd307513f28e3c45021f85c3a4d0c`
- graph_content_hash: `7ee70c556584cf417dda09aa2ff77f50f7b0f7f4e46c2e23846670ac9eb8accc`
- ok_for_seed_mask_gpu: **True**

## Issues

- none

## Per-trait selection (ADR 0012: discovery != G2 input)

| trait | prefilter | discovery | genes | expanded | seed frac | sparsity_ok | strength_cap |
|---|---:|---:|---:|---:|---:|:---:|---:|
| age | 4096 | 1024 | 256 | 9595 | 0.08858780614903596 | False | None |
| sex | 4096 | 44 | 50 | 2356 | 0.01867572156196944 | True | 213.60444111153242 |
| sex_autosome | 4096 | 44 | 50 | 2356 | 0.01867572156196944 | True | 213.60444111153242 |
| tissue | 4096 | 45 | 41 | 2778 | 0.016198704103671708 | True | 126.76232504119935 |

## Overlap

- traits: `['age', 'tissue', 'sex']`
- gene union: 331
- gene pairwise: `{'age_∩_sex': 8, 'age_∩_tissue': 8, 'tissue_∩_sex': 3}`
- CpG union (expanded): 11785
- CpG pairwise: `{'age_∩_sex': 1398, 'age_∩_tissue': 1546, 'tissue_∩_sex': 1186}`
- seed fraction of expanded: `{'age': 0.08858780614903596, 'sex': 0.01867572156196944, 'tissue': 0.016198704103671708}`
- genes with only one seed CpG: `{'age': 35, 'sex': 47, 'tissue': 39}`
- multi-gene CpG count: `{'age': 346, 'sex': 125, 'tissue': 0}`

## G3 matched-random quality

`{'cpg_count_abs_err_max': 265.0, 'cpg_count_abs_err_mean': 2.13595166163142, 'cpg_count_abs_err_median': 0.0, 'cpg_count_abs_err_p90': 0.0, 'fraction_exact_cpg_match': 0.9063444108761329, 'gene_length_bp_used': False, 'gene_role_coverage_used': False, 'n_matched': 331, 'n_seed_genes': 331, 'seed_genes_disjoint_from_matched': True}`
