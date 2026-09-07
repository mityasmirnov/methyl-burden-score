# 7H Phase 4 — reference checkpoints (compact)

ATS cohort (`matrix-hub-age-tissue-sex-full-v1`), split `hub-ats-7e-3fold-v1`,
outer **test** scored once after validation selection. Full fold blobs live
under `artifacts/runs/` (not committed). Raw `summary.json` in this directory
is a large local dump — prefer this file + run metrics.

## P2-G reference refit (`stage0-7h-P2-G-reference`)

Cascade max/max, 15 epochs, gene-linked panel. Product path = `mbs_e2e`.

| Fold | Best ep | Tissue macro-F1 | Age MAE | Sex AUROC |
|-----:|--------:|----------------:|--------:|----------:|
| 0 | 8 | 0.324 | 21.79 | 0.743 |
| 1 | 12 | 0.391 | 13.90 | 0.761 |
| 2 | 5 | 0.422 | 20.99 | 0.737 |
| **mean** | — | **0.379** | **18.89** | **0.747** |

Compare to the prior Stage A `P2-G` screen mean (`mbs_e2e` tissue **0.373**,
age MAE **15.64**). This refit is in the same tissue ballpark; age is weaker
on mean MAE (fold variance is large — do not over-interpret n=3).

## One-hop m-only (matched 16-ep, already in Stage A report)

From `reports/inspection/stage0_7g_gene_only_probe/analysis.md` § matched-budget:

| Arm | Tissue e2e F1 | Age MAE |
|-----|-------------:|--------:|
| `N-light-gene-ablation-m-only-16ep` | **0.372** | **17.7** |
| `N-light-gene-ablation-full-16ep` | 0.351 | 21.8 |

`m_only` remains preferred for one-hop under this budget.

## Next

Phase 3 nine-pack split `hub-nine-pack-3fold-v1` is frozen. **P2-G is not a
product lock** — ATS scalar vs vector was within noise; nine-pack vector RBS
arms are the scale test before any cascade primary. Keep **N-light**
(`m_only`) as a co-equal light MBS finalist for Milestone **12** OOF
(finalists only — do not 5×6 all Stage A arms). Do **not** launch Stage B or
OOF from this ATS checkpoint alone.
