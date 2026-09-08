# Nine-pack scalar vs vector RBS

Generated: `2026-09-07T17:57:28.187880+00:00` (5-combo grid completed `2026-09-08T04:33+02:00`)

Full 5-combo pooling grid now complete at nine-pack scale (34,234 samples), matching
the ATS-scale screen. **P2-G scalar max/max wins outright** — best tissue F1 and best
age MAE of all 5 combos; only scalar mean/max edges it on sex AUROC (+0.022), which is
the lowest-weighted of the three metrics (`sex_loss_weight: 1.0` vs `tissue_loss_weight:
3.0`). Cascade finalist for Milestone 11/12 is **P2-G scalar max/max** — no longer
provisional.

| Arm | folds | Tissue F1 | Age MAE | Sex AUROC |
|-----|------:|----------:|--------:|----------:|
| **P2-G scalar max/max** | 3 | **0.355** | **13.431** | 0.853 |
| scalar mean/max | 3 | 0.330 | 13.496 | 0.875 |
| scalar max/mean | 3 | 0.318 | 17.426 | 0.812 |
| vector mean→max | 3 | 0.335 | 16.448 | 0.780 |
| vector max→max | 3 | 0.333 | 15.408 | 0.834 |
| N-light m-only | 3 | 0.273 | 15.057 | 0.742 |
| N-light m-only, rho_hidden 10→64 | 3 | 0.308 | 14.727 | 0.852 |

Milestone **12** OOF = finalists only (P2-G scalar max/max cascade + N-light).

Full campaign narrative + interpretations: [`analysis.md`](analysis.md).

## Capacity-bottleneck diagnostic (not a cascade substitute)

Widening `rho_hidden_dimension` 10→64 (matching `phi_hidden_dimension`) to test
whether the narrow DeepSet pooling bottleneck explained m-only's widening gap vs.
P2-G at nine-pack scale (vs. near-parity at ATS scale) **helps modestly** vs narrow
m-only (tissue 0.273→0.308, age 15.057→14.727, sex 0.742→0.852) but **does not
close the cascade gap** to P2-G (0.355 / 13.431). Same split/budget/seed/
checkpoint-selection criterion. **Capacity alone is not the bottleneck that
makes one-hop match cascade; keep N-light at `rho_hidden_dimension: 10` for
Milestone 12 OOF** unless a dedicated sweep says otherwise.

## Age MAE: mean vs. median across studies (P2-G baseline)

Aggregate age MAE (13.431) is a mean-of-folds number; it hides how skewed the
per-study distribution is. Across all 190 study×fold rows in the P2-G baseline:

| Statistic | Value |
|---|---|
| mean | 14.16 |
| median | 11.78 |
| p90 | 24.54 |
| min / max | 3.31 / 51.56 |

Median (11.78) sits well below the mean (14.16) — a handful of small-n, hard
studies (colonic mucosa n=3 studies at mean 43.5/median 50.8; brain-DLPFC n=1
at 29.9) pull the mean up. But the median-by-tissue ranking tells the same
story as the mean-by-tissue ranking in the table above: blood/immune tissues
sit at the low end (CD14+ monocyte median 6.8, PBMC 8.2, whole blood 10.5)
and brain regions + colonic mucosa sit at the high end (brain-cerebellum
median 24.9, brain-DLPFC 29.9, colonic mucosa 50.8) regardless of which
statistic you read off. The tissue-conditioning hypothesis doesn't rest on
mean-driven skew — re-check this same mean/median split once the
age-covariates ablation (tissue+sex-conditioned age head) finishes, to see
whether conditioning narrows the *median* gap as well as any mean shift.

## Age-covariates ablation (tissue+sex-conditioned age head): negative result

Full 3-fold run, identical split/budget/seed to the P2-G baseline, only
difference is the age head gets a ground-truth tissue+sex embedding
concatenated in (`model.age_covariates: [tissue, sex]`).

| Metric | P2-G baseline | +tissue/sex conditioning | Δ |
|---|---:|---:|---:|
| Tissue F1 | 0.355 | 0.340 | −0.015 |
| Age MAE | 13.431 | 14.135 | +0.704 (worse) |
| Sex AUROC | 0.853 | 0.834 | −0.019 |

All three metrics moved slightly *against* conditioning, not with it. Per-fold
age MAE: baseline `[12.80, 16.00, 11.49]` vs. conditioned `[12.75, 13.27,
16.38]` — fold 2 flipped from the baseline's best fold to the conditioned
run's worst, which is most of the aggregate swing; folds 0-1 individually
looked promising in isolation, which is exactly why the 3-fold aggregate
(not a single fold) is the number that counts.

Checked whether conditioning at least narrowed the blood-vs-brain MAE gap
(the thing that actually motivated this test) even if the headline MAE
didn't move: blood-like tissues' median MAE 10.2→12.38 (worse), brain-like
tissues' median MAE 20.24→21.67 (worse); gap 10.05→9.29, essentially flat
(~7% narrower, within fold-to-fold noise given the swing above). Overall
median across all study×fold rows: 11.78→11.95 (flat).

**Verdict: does not improve MBS aggregation. Per plan, not adopting
age-covariates conditioning; P2-G scalar max/max (unconditioned) stays the
Milestone 11/12 cascade finalist.** The biological rationale (epigenetic
clocks are tissue-dependent) was sound and the per-tissue heterogeneity is
real, but a linear embedding-concatenation age head evidently isn't the
right lever to exploit it here -- possibly because the shared MBS
representation already implicitly encodes enough tissue signal for the
dense age head to use, or because 3 folds is too few samples of the
effect to detect at this budget. Not planning to pursue a fancier
(nonlinear/FiLM) version of this without new evidence motivating it.

