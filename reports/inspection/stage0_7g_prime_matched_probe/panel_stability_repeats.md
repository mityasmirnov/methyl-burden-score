# Panel stability: `n_repeats=1` (dry-run) vs `n_repeats=5` (production)

Updated: `2026-09-25` — **interim, folds 0–1 of 3** (fold 2 still running).

Compares the committed dry-run panels (`9488b60`, `n_repeats=1`) against the
in-flight production re-run (`n_repeats=5`) for the same folds, same split
(`hub-ats-7e-3fold-v1`), same 65 536-locus universe. Selector unchanged:
`study_grouped_multitask_enet_stability`.

## Result

| Fold | Layer | dry (r=1) | prod (r=5) | shared | Jaccard | dry retained |
|-----:|-------|----------:|-----------:|-------:|--------:|-------------:|
| 0 | `seed_cols`  | 3 410 | 3 409 | 2 939 | **0.757** | 0.862 |
| 0 | `panel_cols` | 32 731 | 33 769 | 31 697 | **0.911** | 0.968 |
| 1 | `seed_cols`  | 3 408 | 3 409 | 2 879 | **0.731** | 0.845 |
| 1 | `panel_cols` | 34 784 | 36 006 | 33 637 | **0.905** | 0.967 |

## Reading it

**Counts are not evidence of convergence — I initially misread them.** The
`n_seed` totals are near-identical across settings (3 410 vs 3 409; 3 408 vs
3 409), which looks like the dry-run had already converged. It had not: only
**~85%** of dry-run seed genes survive at `r=5`, and seed-level Jaccard is
**0.73–0.76**. The similar totals are an artifact of how stability selection
terminates (a frequency threshold against a fixed `max_seeds` budget), so
comparable counts are expected *regardless* of whether the underlying sets
agree. Set overlap had to be computed directly; inferring stability from
counts was wrong.

**The layer that feeds training is much more stable.** `panel_cols` (seed genes
expanded to all their gene-linked CpGs, ADR 0012) retains **~97%** of dry-run
columns at Jaccard **~0.91**. Expansion absorbs most seed-level churn: swapping
a seed gene usually pulls in CpGs that another retained seed already covered.

## Consequence for GATE G1

- Running `n_repeats=5` was **justified** — the seed sets genuinely differ, so
  the dry-run panels were not interchangeable with production ones.
- But the *practical* effect on the G1 gene-set verdict is likely **small**,
  since the trained-on column set is ~97% unchanged. G1's answer ("fold-selected
  panel vs all represented genes") is unlikely to flip on this difference alone.
- **Do not** cite this as G1 evidence yet: fold 2 is still running, and this
  compares panel *composition* only — no model was trained on either panel here,
  so it says nothing about downstream metric impact.
