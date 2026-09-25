# Panel stability: `n_repeats=1` (dry-run) vs `n_repeats=5` (production)

Updated: `2026-09-25` — **CPU chain complete** (panels + classical
`C-mvalue-enetS`, all 3 folds).

Compares the committed dry-run panels (`9488b60`, `n_repeats=1`) against the
production re-run (`n_repeats=5`), same split (`hub-ats-7e-3fold-v1`), same
65 536-locus universe, same selector
(`study_grouped_multitask_enet_stability`).

**Units caveat:** everything below counts **columns** (`seed_cols`,
`panel_cols`), not genes. The "717 genes ∩ all folds" figure quoted elsewhere in
the campaign docs is a *gene*-level count derived after mapping panel columns to
genes, so it is not directly comparable to the column counts here.

## 1. Within-fold reproducibility (r=1 vs r=5, same fold)

| Fold | Layer | dry (r=1) | prod (r=5) | shared | Jaccard | dry retained |
|-----:|-------|----------:|-----------:|-------:|--------:|-------------:|
| 0 | `seed_cols`  | 3 410 | 3 409 | 2 939 | **0.757** | 0.862 |
| 0 | `panel_cols` | 32 731 | 33 769 | 31 697 | **0.911** | 0.968 |
| 1 | `seed_cols`  | 3 408 | 3 409 | 2 879 | **0.731** | 0.845 |
| 1 | `panel_cols` | 34 784 | 36 006 | 33 637 | **0.905** | 0.967 |
| 2 | `seed_cols`  | 3 419 | 3 420 | 2 930 | **0.750** | 0.857 |
| 2 | `panel_cols` | 35 182 | 35 912 | 33 951 | **0.914** | 0.965 |

**Counts are not evidence of convergence — initially misread.** The `n_seed`
totals are near-identical across settings (e.g. 3 410 vs 3 409), which looks
like the dry-run had converged. It had not: only **~85–86%** of dry-run seed
columns survive at r=5 (Jaccard 0.73–0.76). Matching totals are an artifact of
the selector terminating on a frequency threshold against a fixed `max_seeds`
budget, so comparable counts are expected *regardless* of set agreement.

**The layer that feeds training is far more stable:** `panel_cols` (seeds
expanded to all gene-linked CpGs, ADR 0012) retains **~96–97%** at Jaccard
**~0.91**. Expansion absorbs most seed churn — a swapped seed usually pulls in
CpGs another retained seed already covers.

## 2. Cross-fold agreement (do the 3 folds select the same panel?)

| Setting | Layer | ∩ all 3 | ∪ all 3 | ∩/∪ |
|---------|-------|--------:|--------:|-----:|
| dry r=1 | `seed_cols`  | 464 | 7 028 | **0.066** |
| dry r=1 | `panel_cols` | 23 097 | 44 373 | **0.521** |
| prod r=5 | `seed_cols`  | 373 | 7 336 | **0.051** |
| prod r=5 | `panel_cols` | 23 814 | 45 535 | **0.523** |

Two things stand out:

1. **Seed selection is almost entirely fold-specific** — only ~5% of seed
   columns are common to all three folds. Panels overlap much more (~52%), again
   because gene-linked expansion dominates. This is *expected and correct* for
   fold-safe selection (ADR 0011 `internal_fold`: each fold selects on its own
   training data only) — but it means **"the M11 panel" is not one object, it is
   three different panels.**
2. **More repeats did not buy cross-fold agreement.** Going r=1 → r=5 left
   panel-level overlap flat (0.521 → 0.523) and seed-level slightly *lower*
   (0.066 → 0.051). The extra compute bought **within-fold reproducibility**
   (§1), not cross-fold consensus.

## 3. Consequence for GATE G1

- Running `n_repeats=5` was **justified**: the dry-run panels were not
  interchangeable with production ones (~15% of seeds differ).
- But it is **unlikely to flip the G1 verdict**, since the trained-on column set
  is ~96–97% unchanged per fold and cross-fold overlap is essentially identical.
- **New consideration for 10d:** adopting the M11 fold-selected panel as the
  product gene set means shipping a **fold-varying feature set** (~52% common
  across folds). That is fine for a cross-validated *estimate*, but complicates a
  single deployable reference checkpoint — a fixed panel (e.g. the ∩-all-folds
  core, or "all represented genes") sidesteps it. Worth deciding explicitly
  rather than by default.

## 4. Downstream metric impact (classical `C-mvalue-enetS`, r=1 vs r=5)

The classical arm was refit on the r=5 panels, so the composition differences
above can be translated into metrics. Same arm, same folds, same split — the
only change is which panel columns it fit on.

| Metric (3-fold mean) | r=1 (dry) | r=5 (prod) | Δ |
|---|---:|---:|---:|
| age MAE | 8.3095 | 8.3003 | **−0.009** |
| age Pearson r | 0.8941 | 0.8944 | +0.0002 |
| sex AUROC | 0.8905 | 0.8942 | +0.0037 |
| tissue balanced acc | 0.4372 | 0.4395 | +0.0023 |
| tissue macro-F1 | 0.3852 | 0.3883 | +0.0031 |

**`n_repeats=5` did not move the numbers.** Every delta is negligible (age MAE
by 0.1%), and per-fold signs are mixed rather than systematically favouring r=5.
The ~15% seed churn from §1 does not propagate to downstream accuracy, which is
consistent with §1's finding that the trained-on `panel_cols` layer is ~97%
unchanged.

So the documented expectation that "production wants `n_repeats=5`" is now
empirically settled: **r=5 buys within-fold reproducibility, not accuracy.** The
r=1 dry-run was metrically adequate for the G1 decision; the production run
confirms the earlier conclusions rather than revising them. Exploratory panel
work can use r=1 with confidence.

## 5. What this does and does not settle for G1

Settles:
- r=1 vs r=5 is not a live concern (§4).
- The M11 panel is fold-varying, not a single object (§2) — a real input to the
  10d checkpoint decision.

Does **not** settle — G1's actual question is still open:
- G1 asks whether the fold-selected panel **matches or beats "all ~19.6k
  represented genes"** under nested enet. There is **no all-genes comparator
  run**, so no verdict is possible yet. The numbers here characterise the
  fold-selected panel in isolation.
- These are classical `C-mvalue-enetS` numbers, not the neural finalists on the
  same panel, and not `mbs_enet_nested`.
