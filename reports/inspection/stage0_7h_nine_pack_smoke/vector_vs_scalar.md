# Nine-pack scalar vs vector RBS

Updated: `2026-09-08T15:05+02:00`  
Cold 5-combo grid completed `2026-09-08T04:33+02:00`; warm-starts completed `~13:42+02:00`; N-light nested enet 3-fold landed `~14:41+02:00`.

**Architecture-screen readout:** `mbs_e2e` outer **test**, nine-pack split `hub-nine-pack-3fold-v1` (34 234 samples), HM450 only.  
**Product readout (locked):** frozen `rbs_enet` / `mbs_enet(_nested)` co-primary — see §5 and [`milestone-10e-staged-rbs-mbs-training.md`](../../../docs/plans/milestone-10e-staged-rbs-mbs-training.md).  
Full narrative: [`analysis.md`](analysis.md).

---

## Executive summary

| Question | Answer |
|----------|--------|
| Best **cascade** architecture? | **P2-G scalar max/max** (tissue **0.355**, age **13.431**) |
| Does cold vector beat scalar? | **No** — both vector arms trail on tissue + age |
| Was vector just under-trained? | **Partly yes** — LP-FT warm-start closes most of the gap; still no clear win over P2-G |
| Light-model default? | **N-light@64** (0.308 / 14.727 / 0.852) |
| Age + tissue/sex covariates? | **Rejected** (all three metrics worse) |
| Milestone **12** OOF? | **N-light 5×6 first** (queued GPU 2 after S1). Cascade still gated on 10e. **11 does not block** |

**Cascade topology remains P2-G.** Warm vector is an optimization ablation, not a replacement. **Do not 5×6 joint e2e as-is.**

---

## Leaderboard (`mbs_e2e`)

| Arm | folds | Tissue F1 ↑ | Age MAE ↓ | Sex AUROC ↑ |
|-----|------:|------------:|----------:|------------:|
| **P2-G scalar max/max** | 3 | **0.355** | **13.431** | 0.853 |
| vector mean→max **warm** (from scalar mean/max) | 3 | 0.345 | 14.578 | 0.862 |
| vector max→max **warm** (from P2-G) | 3 | 0.342 | **13.406** | 0.851 |
| scalar mean/max | 3 | 0.330 | 13.496 | **0.875** |
| vector mean→max cold | 3 | 0.335 | 16.448 | 0.780 |
| vector max→max cold | 3 | 0.333 | 15.408 | 0.834 |
| scalar max/mean | 3 | 0.318 | 17.426 | 0.812 |
| **N-light m-only (rho=64)** | 3 | 0.308 | 14.727 | 0.852 |
| N-light m-only (rho=10, historical) | 3 | 0.273 | 15.057 | 0.742 |
| P2-G + age covariates | 3 | 0.340 | 14.135 | 0.834 |

---

## 1. Cold 5-combo pooling grid

Matched budget: 15 ep × 3 folds. **P2-G scalar max/max wins** tissue and age; only scalar mean/max edges it on sex (+0.022), the lowest-weighted task (`sex_loss_weight: 1.0` vs `tissue_loss_weight: 3.0`).

| Arm | folds | Tissue F1 | Age MAE | Sex AUROC |
|-----|------:|----------:|--------:|----------:|
| **P2-G scalar max/max** | 3 | **0.355** | **13.431** | 0.853 |
| scalar mean/max | 3 | 0.330 | 13.496 | 0.875 |
| scalar max/mean | 3 | 0.318 | 17.426 | 0.812 |
| vector mean→max | 3 | 0.335 | 16.448 | 0.780 |
| vector max→max | 3 | 0.333 | 15.408 | 0.834 |

**Interpretation.** At nine-pack scale (unlike noisy ATS), scalar vs vector separates cleanly. Cold vector is not primary.

---

## 2. N-light capacity (rho 10→64) — **adopted**

| | rho=10 | rho=64 | Δ |
|--|------:|------:|--:|
| Tissue F1 | 0.273 | **0.308** | +0.035 |
| Age MAE | 15.057 | **14.727** | −0.330 |
| Sex AUROC | 0.742 | **0.852** | +0.110 |

Use **`mbs_e2e`** only. Do **not** cite stale `external_test` 0.194 / 20.6 (under-reports flat tissue F1). Product configs default to `rho_hidden_dimension: 64`. Still trails cascade.

---

## 3. Age-covariates ablation — **rejected**

| Metric | P2-G | + tissue/sex age head | Δ |
|--------|-----:|----------------------:|--:|
| Tissue F1 | 0.355 | 0.340 | −0.015 |
| Age MAE | 13.431 | 14.135 | +0.704 |
| Sex AUROC | 0.853 | 0.834 | −0.019 |

Blood–brain median age gap essentially flat. **Not adopted.**

### Age MAE skew (P2-G baseline, study×fold)

| Statistic | Value |
|-----------|------:|
| mean | 14.16 |
| median | 11.78 |
| p90 | 24.54 |
| min / max | 3.31 / 51.56 |

Blood/immune easy; brain + colonic mucosa hard — under both mean and median.

---

## 4. Vector warm-start (LP-FT) — **both arms done**

Copy scalar encoder → freeze ~4 ep → fine-tune ~11 ep @ 3e-4. Per-fold checkpoints only (no cross-fold leak).

### 4a. max→max warm from P2-G

| Arm | Tissue F1 | Age MAE | Sex AUROC |
|-----|----------:|--------:|----------:|
| P2-G scalar | 0.355 | 13.431 | 0.853 |
| vector max→max cold | 0.333 | 15.408 | 0.834 |
| **vector max→max warm** | **0.342** | **13.406** | **0.851** |

Folds: tissue `[0.348, 0.321, 0.355]`, age `[12.88, 13.18, 14.16]`, sex `[0.949, 0.840, 0.765]`.

### 4b. mean→max warm from scalar mean/max

| Arm | Tissue F1 | Age MAE | Sex AUROC |
|-----|----------:|--------:|----------:|
| scalar mean/max | 0.330 | 13.496 | 0.875 |
| vector mean→max cold | 0.335 | 16.448 | 0.780 |
| **vector mean→max warm** | **0.345** | **14.578** | **0.862** |

Folds age MAE `[12.80, 12.00, 18.94]` — fold 2 volatile (same pattern as age-cov).

**Interpretation.** Warm-start **reliably beats cold vector** and nearly matches P2-G on age/sex for max→max; tissue remains **−0.013 vs P2-G**. Mean→max warm beats its scalar baseline on tissue but loses on age (fold-2 outlier). **Do not switch cascade finalist to vector** — parity-at-best for more complexity. Closes “is vector inherently worse?” → mostly an optimization artifact.

### Optional follow-up (in flight on other GPU)

Dense-gradient stage-1 (`region_pool: mean`) then transplant — `scripts/run_7h_dense_stage1_queue.sh` — because max-pool starves `region_rho` gradients. This is **10e S1 only**, not the full staged recipe and not a 5×6 launch.

---

## 5. Frozen RBS / MBS probes beat `mbs_e2e` — **training lock, not a side note**

Classical probes on region-level RBS (and often gene-pooled MBS) beat jointly trained `mbs_e2e`. This **blocks** launching Milestone 12 on the joint-e2e recipe. Staged RBS→MBS + frozen enet is the pre-OOF plan (10e).

`rbs_linear_probe`: classical (sklearn) probe fit directly on `all_gene_rbs` (`[n_samples, ~15,165 regions]`, pre-gene-pooling), same held-out `test_idx` as `mbs_e2e`, written for every cascade run all along but never compared until now.

| Arm | eval | Tissue F1 | Age MAE | Sex AUROC |
|---|---|---:|---:|---:|
| P2-G scalar max/max | `mbs_e2e` | 0.355 | 13.431 | 0.853 |
| P2-G scalar max/max | `mbs_linear_probe` | 0.347 | **11.567** | 0.816 |
| P2-G scalar max/max | `rbs_linear_probe` | **0.364** | 12.559 | 0.835 |
| vector max/max cold | `mbs_e2e` | 0.333 | 15.408 | 0.834 |
| vector max/max cold | `rbs_linear_probe` | 0.314 | **9.638** | **0.858** |
| vector max/max warm | `mbs_e2e` | 0.342 | 13.406 | 0.851 |
| vector max/max warm | `rbs_linear_probe` | 0.363 | 12.425 | 0.837 |
| scalar mean/max | `mbs_e2e` | 0.330 | 13.496 | 0.875 |
| scalar mean/max | `rbs_linear_probe` | **0.368** | 12.059 | 0.859 |

Consistent across every arm checked: `rbs_linear_probe` beats `mbs_e2e` on age MAE every time (0.9-5.8 years), and beats or ties it on tissue F1 in 3/4 arms. The best tissue F1 (0.368) and best age MAE (9.638) of the whole campaign both come from `rbs_linear_probe`, not any `mbs_e2e` arm. `mbs_linear_probe` (classical on gene-pooled MBS, not region-level RBS) also beats `mbs_e2e` on age MAE — so part of the gap is "classical regression beats a jointly-trained multi-task neural head," part is "raw region-level features carry more signal than gene-pooled ones," stacking in the RBS case.

---

## Outlook

1. **Topology** locked for the e2e screen: P2-G cascade + N-light@64 encoder. Rank under **nested enet**.  
2. **One-hop smokes done** (G0 ≫ G1; seeds 42/43/44 distinct).  
3. **N-light nested enet 3-fold 0.368 / 9.88 / 0.803** vs e2e 0.308 / 14.73 / 0.852.  
4. **Milestone 12 N-light 5×6** queued on GPU 2 after dense S1 (skip transplants).  
5. Cascade 5×6 still waits on 10e S1–S4. Extra traits with **n≥600** = freeze-reuse. Encoder stays gene-invariant for later EPIC/ONT.
