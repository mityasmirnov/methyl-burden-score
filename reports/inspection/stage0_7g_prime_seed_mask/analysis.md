# 7G′ age-primary seed-mask screen

Selection: validation age MAE primary; tissue F1 secondary; sex AUROC tertiary.
Seed source: **internal_fold** (ADR 0011). Discovery CpGs ≠ G2 input (ADR 0012).
P2-G topology is a reference, not a lock.

Folds: [0]; seeds: [42, 43]; K=256.

Arms: G0 all-gene control; G1 head masks; G2 expanded gene CpGs+masks;
G3 matched random; C0 classical all-gene; C2 classical on G2 expanded CpGs.

graph_content_hash: `7ee70c556584cf417dda09aa2ff77f50f7b0f7f4e46c2e23846670ac9eb8accc`
panel_hash: `ef6cd307513f28e3c45021f85c3a4d0c`
configured_traits: `[{'autosome_control': False, 'id': 'age', 'role': 'primary'}, {'autosome_control': False, 'id': 'tissue', 'role': 'secondary'}, {'autosome_control': True, 'id': 'sex', 'role': 'auxiliary'}]`

## Per-trait discovery vs expanded (fold 0)

| trait | prefilter | discovery CpGs | seed genes | unique expanded | edges | seed frac | sparsity_ok |
|---|---:|---:|---:|---:|---:|---:|:---:|
| age | 4096 | 1024 | 256 | 9595 | 9941 | 0.08858780614903596 | False |
| sex | 4096 | 44 | 50 | 2356 | 2481 | 0.01867572156196944 | True |
| sex_autosome | 4096 | 44 | 50 | 2356 | 2481 | 0.01867572156196944 | True |
| tissue | 4096 | 45 | 41 | 2778 | 2778 | 0.016198704103671708 | True |

## Overlap (configured traits)

- traits: `['age', 'tissue', 'sex']`
- gene set sizes: `{'age': 256, 'sex': 50, 'tissue': 41}`
- gene union: 331
- gene pairwise: `{'age_∩_sex': 8, 'age_∩_tissue': 8, 'tissue_∩_sex': 3}`
- expanded CpG set sizes: `{'age': 9595, 'sex': 2356, 'tissue': 2778}`
- expanded CpG union: 11785
- CpG pairwise: `{'age_∩_sex': 1398, 'age_∩_tissue': 1546, 'tissue_∩_sex': 1186}`
- seed fraction of expanded: `{'age': 0.08858780614903596, 'sex': 0.01867572156196944, 'tissue': 0.016198704103671708}`
- gene-role coverage: `{'age': {'five_prime': 389, 'gene_body': 5457, 'promoter_core': 1855, 'promoter_proximal': 1773, 'three_prime': 467}, 'sex': {'five_prime': 36, 'gene_body': 1647, 'promoter_core': 398, 'promoter_proximal': 301, 'three_prime': 99}, 'tissue': {'five_prime': 49, 'gene_body': 1990, 'promoter_core': 307, 'promoter_proximal': 319, 'three_prime': 113}}`
- genes with only one discovery CpG: `{'age': 35, 'sex': 47, 'tissue': 39}`
- multi-gene CpG count: `{'age': 346, 'sex': 125, 'tissue': 0}`

## G3 matched-random quality

`{'cpg_count_abs_err_max': 265.0, 'cpg_count_abs_err_mean': 2.13595166163142, 'cpg_count_abs_err_median': 0.0, 'cpg_count_abs_err_p90': 0.0, 'fraction_exact_cpg_match': 0.9063444108761329, 'gene_length_bp_used': False, 'gene_role_coverage_used': False, 'n_matched': 331, 'n_seed_genes': 331, 'seed_genes_disjoint_from_matched': True}`

See `/data/projects/methyl-burden-score/reports/inspection/stage0_7g_prime_seed_mask/summary.json` and `panel_audit.md`.
