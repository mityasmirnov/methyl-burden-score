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
| N-light m-only (rho=10, historical) | 3 | 0.273 | 15.057 | 0.742 |
| **N-light m-only (rho=64, default)** | 3 | **0.308** | **14.727** | **0.852** |

All rows are **`mbs_e2e`** outer test (not `external_test`, which under-reports tissue F1
for flat runs and previously caused a false “rho64 worse” reading of 0.194 / 20.6).

Milestone **12** OOF = finalists only (P2-G scalar max/max cascade + **N-light@64**).

Full campaign narrative + interpretations: [`analysis.md`](analysis.md).

## Capacity: rho 10→64 adopted as N-light default

Widening `rho_hidden_dimension` 10→64 (matching `phi_hidden_dimension`), evaluated
via primary **`mbs_e2e`** on `stage0-7h-nine-pack-m-only-wide-f*`:

| | rho=10 | rho=64 | Δ |
|--|------:|------:|--:|
| Tissue F1 | 0.273 | **0.308** | +0.035 |
| Age MAE | 15.057 | **14.727** | −0.330 |
| Sex AUROC | 0.742 | **0.852** | +0.110 |

**Adopted as the N-light default** (2026-09-08) in product configs + `flat_region`
code fallbacks. Still trails P2-G cascade (0.355 / 13.431) — capacity helps one-hop
but does not replace cascade. Do **not** cite the stale `external_test` 0.194/20.6
figures for this comparison.

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

## Vector (region_hidden) warm-started from P2-G scalar checkpoint

Motivation: from-scratch vector max/max scored 0.333 F1 / 15.408 MAE / 0.834
AUROC vs. P2-G scalar's 0.355 / 13.431 / 0.853 -- worse on every metric. But
`gene_rho` (vector's only extra learned module) was trained jointly from
random init alongside the CpG/region encoder in that run, which per LP-FT
(Kumar et al. 2022) risks the random head's early noisy gradients distorting
an otherwise-good encoder. New mechanism built in `cascade_loop.py`
(`_load_encoder_warm_start`, `warm_start_encoder_checkpoint` /
`freeze_encoder_epochs` / `fine_tune_learning_rate` config fields, unit-tested
in `test_stage0_7g_tissue_probe.py::test_cascade_warm_start_*`): copy the
converged P2-G `cpg_encoder`/`region_type_embedding`/`region_encoder`/
`region_rho` into a fresh `region_hidden` model (per-fold, from that same
fold's own P2-G checkpoint -- never a different fold, which would leak that
fold's train split into this fold's initialization), freeze it for 4 epochs
while `gene_rho` + trait heads linear-probe, then unfreeze and fine-tune
everything jointly at 3e-4 for the remaining 11 epochs (15 total, matching
every other arm's budget).

**Full 3-fold result: warm-start closes almost the entire scalar-vs-vector gap.**

| Arm | folds | Tissue F1 | Age MAE | Sex AUROC |
|-----|------:|----------:|--------:|----------:|
| P2-G scalar max/max (baseline) | 3 | 0.355 | 13.431 | 0.853 |
| vector max/max, **from scratch** | 3 | 0.333 | 15.408 | 0.834 |
| vector max/max, **warm-started** | 3 | **0.342** | **13.406** | **0.851** |

Per-fold: tissue F1 `[0.348, 0.321, 0.355]`, age MAE `[12.88, 13.18, 14.16]`,
sex AUROC `[0.949, 0.840, 0.765]`. Warm-started vector is essentially at
parity with P2-G scalar (age MAE even marginally better, 13.406 vs 13.431;
tissue F1 0.013 lower; sex AUROC 0.002 lower — all within fold-to-fold
noise) and clearly beats from-scratch vector on every metric (tissue F1
+0.009, age MAE −2.0 years, sex AUROC +0.017). This strongly supports the
LP-FT hypothesis: vector's earlier underperformance was substantially an
optimization artifact of training a random-init `gene_rho` jointly with the
encoder from scratch, not a fundamental limitation of gene-level vector
aggregation vs. scalar pooling. **Not adopting vector as the new finalist**
over P2-G scalar (no clear win, added complexity for parity at best), but
this closes the "is vector inherently worse" question raised by the
original from-scratch grid — it isn't; it just needed a better training
recipe. P2-G scalar max/max remains the Milestone 11/12 cascade finalist.

**Vector mean/max warm-started from scalar mean/max** (second warm-start
arm, same LP-FT recipe): tissue F1 **0.345**, age MAE **14.578**, sex AUROC
**0.862**, vs. scalar mean/max baseline (0.330/13.496/0.875) and
from-scratch vector mean/max (0.335/16.448/0.780). Beats from-scratch
vector mean/max clearly on every metric (age MAE −1.87 years, sex AUROC
+0.082). Vs. the scalar mean/max baseline: tissue F1 better (+0.015), age
MAE worse (+1.08) but driven by one volatile fold (per-fold age MAE
`[12.80, 12.00, 18.94]` — fold 2 is an outlier, same pattern seen in the
age-covariates ablation's fold 2), sex AUROC essentially flat (−0.013).
Net: same conclusion as max/max — warm-start reliably closes most or all
of the from-scratch gap, confirming this isn't a max/max-specific fluke,
but doesn't clearly surpass scalar. **Both warm-start arms corroborate:
vector was underrated by the original from-scratch grid, but scalar
remains the safer, no-added-complexity finalist choice.**

### Known limitation of this warm-start (not yet addressed)

Even P2-G's own training never gave `region_rho` a *dense* gradient signal:
`scalar_rbs`'s "MBS" is `max_pool(RBS scores in gene)`, and max-pooling
routes gradient only to the single argmax region per gene per sample --
every other region gets exactly zero gradient that step. So the checkpoint
being warm-started from was itself trained under a gradient-starved regime
for `region_rho`, the whole 15 epochs. A more thorough fix would be a
dedicated stage-1 pretraining pass with `region_pool: mean` (dense gradient
to every region) purely to train the encoder well, then transplant *that*
into the max-pooled stage-2 config (max-pooling is still likely the right
choice for the *final* aggregation -- burden scores want "is there at least
one impaired region," which mean-pooling would dilute; `scalar_max_mean`
(region_pool: mean) already scored worse than `scalar_max_max`, 0.318 vs
0.355, consistent with this). Not yet built -- proposed as a follow-up arm
if the simpler warm-start above doesn't close the gap.

## RBS-only classical probe vs. MBS-based prediction (existing eval, newly surfaced)

The cascade pipeline already writes `rbs_linear_probe` for every arm: a
classical (sklearn) probe fit directly on `all_gene_rbs`
(`[n_samples, n_regions]`, ~15,165 columns -- every region's RBS score,
*before* any gene-level pooling), trained on the same `train_idx` and
scored on the same held-out `test_idx` as `mbs_e2e`. This has been computed
all along but not compared against the neural end-to-end head until now.

| Arm | eval | Tissue F1 | Age MAE | Sex AUROC |
|---|---|---:|---:|---:|
| P2-G scalar max/max | `mbs_e2e` (neural, joint) | 0.355 | 13.431 | 0.853 |
| P2-G scalar max/max | `mbs_linear_probe` (classical, MBS) | 0.347 | **11.567** | 0.816 |
| P2-G scalar max/max | `rbs_linear_probe` (classical, RBS) | **0.364** | 12.559 | 0.835 |
| vector max/max, from scratch | `mbs_e2e` | 0.333 | 15.408 | 0.834 |
| vector max/max, from scratch | `rbs_linear_probe` | 0.314 | **9.638** | **0.858** |
| vector max/max, warm-started | `mbs_e2e` | 0.342 | 13.406 | 0.851 |
| vector max/max, warm-started | `rbs_linear_probe` | 0.363 | 12.425 | 0.837 |
| scalar mean/max | `mbs_e2e` | 0.330 | 13.496 | 0.875 |
| scalar mean/max | `rbs_linear_probe` | **0.368** | 12.059 | 0.859 |

**Consistent pattern across every arm tested: a classical probe on raw
region-level RBS beats the neural `mbs_e2e` head on age MAE, every single
time** (by 0.9-5.8 years), and beats or ties it on tissue F1 in 3 of 4
arms. The best tissue F1 in this entire campaign (0.368) and the best age
MAE (9.638, from-scratch vector) both come from `rbs_linear_probe`, not
from any `mbs_e2e` arm we've been optimizing. `mbs_linear_probe` (classical
probe on the *gene-pooled* MBS, not RBS) also beats `mbs_e2e` on age MAE
in every arm checked -- so part of this gap is "classical regression beats
a jointly-trained multi-task neural head," and part is specifically "raw
region-level features carry more signal than gene-pooled ones," stacking
in the RBS case.

**Implication:** the neural `mbs_e2e` head has been our primary evaluation
and finalist-selection metric this whole campaign, but it is not the best
predictor available from these architectures' own features. This directly
supports investing more in the region-level (RBS) pathway rather than
treating gene pooling as strictly necessary for the best score -- consistent
with the dense-gradient `region_rho` limitation noted above (a better-
trained RBS encoder plausibly widens this gap further). **Not yet decided:**
whether Milestone 12 OOF should report/deploy `rbs_linear_probe` (or an
RBS+MBS fusion) as a co-primary readout alongside `mbs_e2e`, given it's
already computed for free. Flagging for discussion, not changing the
primary evaluation unilaterally.

