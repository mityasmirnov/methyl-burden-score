# Plan: Matched 16-epoch promotion screen (7G′ Stage A)

Status: **in progress** — current GPU gate. Age-primary seed-mask stays
**blocked** until the decision rules below fire.

Parent: [`milestone-7g-prime-stage-a-deeprvat-screen.md`](milestone-7g-prime-stage-a-deeprvat-screen.md).
Report: `reports/inspection/stage0_7g_gene_only_probe/analysis.md`.

## Goal

Match epoch budgets before declaring scalar cascade better than vector or
one-hop. Keep **P2-G as the current reference, not a pooling lock.** Fixed
`rbs_enet` (α=0.1, l1_ratio=0.5, no scaler) is diagnostic only — age collapse
with stable sex on vector RBS is a sparse-penalty failure, not evidence the
vector representation is weak.

## Locked decisions

| Choice | Decision |
|--------|----------|
| Ceiling | 16 epochs; patience 5; tissue-primary checkpoint |
| Loss | Keep tissue-primary (`λ_t=3`, `λ_a=0.3`, `λ_s=1`) |
| Panel | `explicit_only` 51 375 CpGs |
| Microbatch | one-hop batch 256 / 16M tokens; cascade batch **128** |
| References | Keep P2-G / P4-G (15 ep); do not rerun |
| Vector max→max | Do **not** promote yet |
| Nested enet | Post-hoc CPU; train-fold StandardScaler + inner-val α/l1 |

## Queue (14 GPU-fold jobs on GPU 0)

1. `N-light-gene-max` folds 1–2 (fold 0 already selected epoch 16)
2. `N-light-gene-mean` all 3 folds (new `…-16ep` prefix)
3. Scalar `mean→max` and `max→mean` all 3 folds (new `…-16ep` run IDs)
4. Vector `mean→max` all 3 folds (new `…-16ep` run ID)

Driver: `scripts/run_7g_16ep_promotion.sh` (resume: `…_resume.sh`).
**Nested / fixed elastic-net is post-hoc only** — do not run
`eval_mbs_enet_from_scores.py` inside the GPU queue; refresh the report
from e2e + linear/Ridge probes first, then enet offline.

## Decision rules

- If one-hop max ≈ P2-G across three folds → preferred smaller DeepRVAT-like architecture.
- If vector mean→max improves age/sex e2e but still loses after gene pooling → typed-RBS aggregation, not scalar MBS.
- If scalar mixed pooling closes the gap to P2 → retain the full 2×2 pooling result.
- If nothing beats P2-G or raw-CpG baselines → stop architecture sweeps; start age-primary seed-gene.

## Non-goals

- Launching seed-mask GPU training before the rules fire
- Stage B CpG-panel GPU / Milestone 7
- Rerunning P2/P4 for exact 16-ep equality
- Promoting vector max→max
