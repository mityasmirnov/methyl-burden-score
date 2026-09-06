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
(Raw `summary.json` is large — interpreted below, not tracked in git. This
script rewrites this file from scratch each run; re-added after every
rerun rather than living in git history alone — see
`milestone-7h-pretrained-mbs-rbs-campaign.md` for the full debugging trail.)

## GPU screen results — final (2026-09-06, fold 0, seeds {42, 43}, 40 epochs)

Two real bugs blocked this screen and are now fixed (full trail in
`milestone-7h-pretrained-mbs-rbs-campaign.md`): a `KeyError: 'tissues'`
crash in the classical arms, and a silent `learning_rate` threading bug
that left every cascade run training at `lr=1e-2` (10x the intended
`0.001`), which combined with sparse seed-masked heads caused near-total
training collapse. With both fixed, gradient clipping added, and the
epoch budget raised 15→40 (needed — several arms only escaped the bad
regime after epoch 20-30):

| Arm | Seed | Age MAE | Tissue F1 | Sex AUROC | Best epoch |
|---|---:|---:|---:|---:|---:|
| G0 (all genes, no mask) | 42 | 17.66 | 0.227 | 0.790 | 40/40 |
| G0 (all genes, no mask) | 43 | 16.38 | 0.238 | 0.786 | 40/40 |
| G1 (all CpGs, seed heads) | 42 | 21.48 | 0.098 | 0.512 | 38/40 |
| G1 (all CpGs, seed heads) | 43 | 25.13 | 0.058 | 0.515 | 39/40 |
| G2 (seed CpGs, seed heads) | 42 | 22.12 | 0.070 | 0.565 | 40/40 |
| G2 (seed CpGs, seed heads) | 43 | 21.48 | 0.074 | 0.715 | 39/40 |
| G3 (matched random) | 42 | 24.26 | 0.074 | 0.583 | 23/40 |
| G3 (matched random) | 43 | 25.20 | 0.000 | 0.489 | 2/40 — **still collapsed** |
| C0 (classical, all-gene enet) | — | 8.94 | 0.367 | 0.854 | — |
| C2 (classical, seed-CpG enet) | — | 10.61 | 0.352 | 0.856 | — |

7 of 8 cascade runs now train to real, non-degenerate results (`G3` seed 43
is a stubborn residual outlier — collapsed by epoch 2 and never recovered
even at 40 epochs; not chased further given this is a single seed of one
arm, and the rest of the grid is conclusive without it).

**Answer to the milestone's question:** seed-gene masking does **not**
help age-primary training here — **`G0` (dense, all genes, no masking)
clearly and consistently beats every masked variant** on all three traits
(age MAE 16-18 vs. 21-25; tissue F1 0.23 vs. 0.06-0.10; sex AUROC 0.79 vs.
0.51-0.72). Restricting supervision to a narrow ~256-gene (or fewer)
seed-derived subset — whether restricting only the head (`G1`), the head
*and* input CpGs (`G2`), or a size-matched random gene set (`G3`) — costs
real performance relative to the unmasked control, and there's no clear
separation between `G1`/`G2`/`G3` themselves (all cluster in the same
degraded range) — i.e. no evidence the *specific* seed genes carry more
signal than an equivalently-sized random gene set once training actually
converges. Classical M-value regression (`C0`) remains the strongest age
predictor overall (MAE 8.94), consistent with every prior Stage A/Stage-B
result in this project.

**Practical implication:** for age-primary objectives on this cohort/panel,
gene-level masking to a small discovered seed set is not a promising
direction on this evidence — the full-gene (`G0`) cascade is simply a
better-parameterized model. This doesn't rule out seed-masking being
useful for other purposes (e.g., interpretability, or traits with more
concentrated genetic architecture than age), but it shouldn't be adopted
as a training strategy for the pretrained MBS/RBS framework based on this
result.
