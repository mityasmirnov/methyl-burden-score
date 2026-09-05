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
curves — interpreted below, not tracked in git. Note: `run_7g_prime_seed_mask.py`
rewrites this file from scratch each run, so this section is re-added after
every rerun rather than living only in git history.)

## GPU screen results (2026-09-05, fold 0, seeds {42, 43})

First launch (weekend supervisor) crashed on `C0` seed 42 with
`KeyError: 'tissues'` in `classical_mvalue.fit_eval_mvalue_fold` — fixed by
threading `class_names` through `_phenotype_arrays()` in
`run_7g_prime_seed_mask.py` / `run_7g_prime_stage_b.py` (`d00a7b1`).

**Classical arms (unaffected by any of this, both runs identical):**

| Arm | Age MAE | Age R² | Tissue F1 | Sex AUROC |
|---|---:|---:|---:|---:|
| `C0` (all-gene enet) | 8.94 | 0.779 | 0.367 | 0.854 |
| `C2` (seed-CpG enet) | 10.61 | 0.704 | 0.352 | 0.856 |

**Cascade arms, run 1 (no gradient clipping) — 7 of 8 collapsed:** age
MAE≈25 (≈ predicting the training mean), tissue F1≈0, sex AUROC=0.5, flat
across all 15 epochs; `encoder_grad_norm` went from normal (0.03–0.6) at
epoch 1 to ~1e-9–1e-13 by epoch 3 in every collapsed run, including `G0`
(no seed-masking at all) — pointing to a generic cascade-trainer
instability, not a `SeedMaskedLinearHead` bug specifically. Root cause:
`train_cascade_on_arrays` had no gradient clipping at all (unlike the
flat-baseline trainer, which defaults to `gradient_clip_norm=2.0`), combined
with this screen's age-primary loss weights (`age_loss_weight=1.0` vs. the
`0.3` used by every previously-validated tissue-primary cascade arm).

**Fix applied:** added `gradient_clip_norm` (default 2.0, matching
`loop.py`) to `train_cascade_on_arrays`, threaded through both call sites.
All 9 existing cascade unit tests still pass.

**Cascade arms, run 2 (gradient clipping added) — `G0` fixed, `G1`/`G2`/`G3`
still collapse:**

| Arm | Seed | Age MAE | Tissue F1 | Sex AUROC | Status |
|---|---:|---:|---:|---:|---|
| G0 (all genes, no mask) | 42 | 13.07 | 0.142 | 0.704 | **fixed** |
| G0 (all genes, no mask) | 43 | 13.99 | 0.170 | 0.822 | **fixed** |
| G1 (all CpGs, seed heads) | 42 | 25.08 | 0.0001 | 0.500 | still collapsed |
| G1 (all CpGs, seed heads) | 43 | 25.10 | 0.0001 | 0.500 | still collapsed |
| G2 (seed CpGs, seed heads) | 42 | 25.08 | 0.0001 | 0.500 | still collapsed |
| G2 (seed CpGs, seed heads) | 43 | 25.10 | 0.0001 | 0.500 | still collapsed |
| G3 (matched random) | 42 | 25.09 | 0.0000 | 0.500 | still collapsed |
| G3 (matched random) | 43 | 25.08 | 0.0034 | 0.500 | still collapsed |

Gradient clipping fully fixed `G0` (both seeds now train normally, matching
classical-arm ballpark). **All 6 seed-masked runs (`G1`/`G2`/`G3` × 2 seeds)
still collapse identically to before** — same flat epoch-1-through-15
trajectory, same `encoder_grad_norm` decay to ~0 by epoch 3-6. This isolates
the remaining failure specifically to `SeedMaskedLinearHead`
(`src/mbs/models.py:611`), not the generic cascade-training instability
clipping just fixed.

Checked and ruled out: the seed masks themselves are not degenerate (`age`
256/2646 genes active, `sex` 50/2646, `tissue` 41/2646 per class row across
47 classes for G1/G2; `G3`'s matched-random masks are similarly sized) — so
this isn't a masks-are-all-zero bug.

**Working (unconfirmed) hypothesis:** `SeedMaskedLinearHead.gene_weight` is
zero-initialized (`nn.Parameter(torch.zeros(n_outputs, n_genes))`), and the
effective weight is `gene_weight * seed_mask`. The gradient into the shared
encoder through this head is `d(loss)/d(mbs) = d(loss)/d(output) @ weight.T`,
which is exactly zero at initialization since `weight = 0`. Combined with
AdamW `weight_decay` acting on `gene_weight` every step, a head whose only
active parameters are a small masked subset (41–331 of 2646 genes) may have
too weak/noisy a gradient signal to escape zero net drift, while `G0`'s
dense 2646-gene head has a much stronger aggregate gradient that overcomes
decay easily. Not yet confirmed — would need e.g. a small-scale random init
for `gene_weight`, or excluding masked-head params from weight decay, then
a targeted rerun on `G1` alone before spending a full 8-run grid again.

**Recommendation:** do not treat `G1`/`G2`/`G3` results as valid — the
G0-vs-seed-masking comparison this milestone exists to answer is still
unanswered. `G0`'s numbers (now trained normally) and `C0`/`C2` are usable.
Next step is architectural (fix `SeedMaskedLinearHead` init/optimizer
interaction), not another blind rerun.
