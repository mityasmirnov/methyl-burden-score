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
(Raw `summary.json` is 218 MB, mostly per-study confusion matrices and ROC
curves — interpreted below, not tracked in git.)

## GPU screen results (2026-09-05, fold 0, seeds {42, 43})

First launch (weekend supervisor) crashed on `C0` seed 42 with
`KeyError: 'tissues'` in `classical_mvalue.fit_eval_mvalue_fold`
(`_phenotype_arrays()` in `run_7g_prime_seed_mask.py` / `run_7g_prime_stage_b.py`
never populated a string-name `ph["tissues"]` array). Fixed by threading
`class_names` through both scripts; relaunched the full grid and it completed
cleanly end to end.

**Classical arms trained normally:**

| Arm | Age MAE | Age R² | Tissue F1 | Sex AUROC |
|---|---:|---:|---:|---:|
| `C0` (all-gene enet) | 8.94 | 0.779 | 0.367 | 0.854 |
| `C2` (seed-CpG enet) | 10.61 | 0.704 | 0.352 | 0.856 |

**Cascade arms (G0–G3, `mbs_e2e`) mostly did not train — do not use for a
G0-vs-G1/G2/G3 decision yet:**

| Arm | Seed | Age MAE | Age R² | Tissue F1 | Sex AUROC | Status |
|---|---:|---:|---:|---:|---:|---|
| G0 (all genes, no mask) | 42 | 25.10 | −0.03 | 0.000 | 0.500 | **collapsed** |
| G0 (all genes, no mask) | 43 | 13.93 | 0.474 | 0.304 | 0.720 | trained |
| G1 (all CpGs, seed heads) | 42 | 25.08 | ~0 | 0.000 | 0.500 | **collapsed** |
| G1 (all CpGs, seed heads) | 43 | 25.10 | ~0 | 0.000 | 0.500 | **collapsed** |
| G2 (seed CpGs, seed heads) | 42 | 25.08 | ~0 | 0.000 | 0.500 | **collapsed** |
| G2 (seed CpGs, seed heads) | 43 | 25.10 | ~0 | 0.000 | 0.500 | **collapsed** |
| G3 (matched random) | 42 | 25.09 | ~0 | 0.000 | 0.500 | **collapsed** |
| G3 (matched random) | 43 | 25.08 | ~0 | 0.000 | 0.500 | **collapsed** |

**7 of 8 cascade runs collapsed** to the random/constant-prediction baseline
(tissue F1 ≈ 0, sex AUROC = 0.5, age MAE ≈ 25 ≈ predicting the training
mean) and never recovered across all 15 epochs. `val_history` shows this is
not a checkpoint-selection artifact — age MAE/tissue F1/sex AUROC are flat
at the random-baseline level from epoch 1 through 15 in every collapsed run.
Per-batch logs show `encoder_grad_norm` starting at a normal order of
magnitude (0.03–0.6) at epoch 1 but collapsing to ~1e-9–1e-13 by epoch 3 in
every collapsed run — consistent with an early large update saturating the
encoder (dead-gradient region) rather than random noise. Collapse happens
**even for `G0` (no seed-masking at all)**, so this is not specific to
`SeedMaskedLinearHead` — it reproduces with plain cascade training under
this config's loss weights.

**Likely cause:** `configs/experiment/stage0_7g_prime_seed_mask.yaml` sets
`age_loss_weight: 1.0` (age-primary, vs. `lambda_age: 0.3` on every
previously-validated tissue-primary cascade arm) with **no
`gradient_clip_norm`** and **no age-target standardization** — neither
exists in `train_cascade_on_arrays` (`src/mbs/training/cascade_loop.py`) at
all, unlike the flat-baseline trainer (`loop.py`) which has both. Raw
unstandardized age (years) driving the dominant loss term with unclipped
gradients is a plausible mechanism for the saturating first-few-epochs
blowup. This is the **first age-primary GPU run** in the project, so no
previously-committed result depends on this loss-weight regime — earlier
cascade arms (P2-G, N-cascade-*) are tissue-primary (`lambda_tissue: 3.0`)
and did not hit this failure mode.

**Recommendation:** do not treat this run as the seed-mask decision. Add
gradient clipping (and ideally train-fold age standardization, matching
`loop.py`'s `target_standardization: train_fold`) to
`train_cascade_on_arrays`, then rerun G0–G3 before drawing any conclusion
about seed-masking's effect on age-primary training.
