# Milestone 10 — nine-pack campaign analysis

Updated: `2026-09-08T13:07+02:00`

- Matrix: `matrix-hub-nine-pack-virtual-v1` (virtual multi-store)
- Split: `hub-nine-pack-3fold-v1` (**34 234** samples)
- Platform: **HM450 only** — no cross-platform claim
- Primary readout: **`mbs_e2e`** on outer **test** (means ± population SD across 3 folds unless noted)

**How to read this file:** newest / decision-blocking results first; each section ends with an explicit **Interpretation**. Machine snapshots: [`summary.json`](summary.json), [`vector_vs_scalar.md`](vector_vs_scalar.md).

---

## Live status (GPU 0)

| Job | Status | Detail |
|-----|--------|--------|
| **Vector mean→max warm-start** | **running** | Fold 0 done; fold 1 ~epoch 12/15 (LP-FT from scalar mean/max ckpt) |
| Vector max→max warm-start | **done** (3/3) | See §4 — improves vs cold; still trails P2-G tissue |
| Age-covariates ablation | done (rejected) | §2 |
| One-hop correctness smokes | **queued** | Waiter blocked until warm-start queue exits; import fix landed |
| Keeper soft-stop | idle | Prior vector→onehop→B.4 keeper finished; ATS B.4 already done |

---

## Executive decisions (current)

| Decision | Status | Basis |
|----------|--------|-------|
| Cascade finalist | **P2-G scalar max/max** | Best tissue F1 **and** best age MAE in cold 5-combo grid |
| Light-model finalist | **N-light m-only (`rho_hidden=64`)** | Wide `mbs_e2e` beat rho=10 on all three metrics |
| Cold vector RBS | **Does not beat P2-G** | Both cold vector arms trail tissue + age |
| Warm vector max→max | **Closer, still not primary** | 0.342 / **13.41** / 0.851 vs P2-G 0.355 / 13.43 / 0.853 |
| Age head + tissue/sex covariates | **Not adopted** | All three headlines slightly worse |
| Milestone **12** OOF | **Finalists only** | P2-G + N-light@64 — not full arm matrix / not all 9 packs |

---

## Leaderboard (`mbs_e2e`, nine-pack test)

| Arm | folds | Tissue F1 ↑ | Age MAE ↓ | Sex AUROC ↑ | Note |
|-----|------:|------------:|----------:|------------:|------|
| **P2-G scalar max/max** | 3 | **0.355** | **13.431** | 0.853 | **Cascade finalist** |
| vector max→max **warm** | 3 | 0.342 | **13.406** | 0.851 | Closes most of age gap; −0.013 tissue vs P2-G |
| scalar mean/max | 3 | 0.330 | 13.496 | **0.875** | Best sex only |
| vector mean→max cold | 3 | 0.335 | 16.448 | 0.780 | |
| vector max→max cold | 3 | 0.333 | 15.408 | 0.834 | |
| scalar max/mean | 3 | 0.318 | 17.426 | 0.812 | |
| **N-light rho=64** | 3 | 0.308 | 14.727 | 0.852 | **Light finalist** |
| N-light rho=10 | 3 | 0.273 | 15.057 | 0.742 | historical |
| P2-G + age covariates | 3 | 0.340 | 14.135 | 0.834 | rejected |
| vector mean→max **warm** | 1+ | 0.364* | 12.80* | 0.941* | *fold 0 only — in flight |

**Interpretation.** At nine-pack scale the story is stable: **scalar P2-G remains the cascade pick**. Warm-starting vector max→max from the P2-G checkpoint **mostly fixes age** (15.4→13.4, matching P2-G) and lifts tissue (0.333→0.342) but **does not overtake** P2-G’s 0.355 tissue. Treat warm vector as an optimization ablation, not a new finalist, unless mean→max warm 3-fold clearly beats P2-G (fold 0 alone is not enough).

---

## 1. Scalar vs vector RBS (cold) — scale screen

Matched budget: **15 epochs × 3 folds**.

| Arm | folds | Tissue F1 | Age MAE ↓ | Sex AUROC |
|-----|------:|----------:|----------:|----------:|
| **P2-G scalar max/max** | 3 | **0.355** (±0.044) | **13.431** (±1.894) | 0.853 (±0.071) |
| scalar mean/max | 3 | 0.330 (±0.013) | 13.496 (±1.682) | **0.875** (±0.072) |
| scalar max/mean | 3 | 0.318 (±0.031) | 17.426 (±1.797) | 0.812 (±0.085) |
| vector mean→max | 3 | 0.335 (±0.033) | 16.448 (±1.724) | 0.780 (±0.056) |
| vector max→max | 3 | 0.333 (±0.039) | 15.408 (±3.806) | 0.834 (±0.050) |

Fold detail (tissue / age / sex):

| Arm | fold 0 | fold 1 | fold 2 |
|-----|--------|--------|--------|
| P2-G | 0.364 / 12.80 / 0.950 | 0.297 / 16.00 / 0.830 | 0.404 / 11.49 / 0.780 |
| scalar mean/max | 0.331 / 13.84 / 0.927 | 0.314 / 11.29 / 0.924 | 0.346 / 15.37 / 0.774 |
| scalar max/mean | 0.310 / 15.40 / 0.906 | 0.284 / 17.12 / 0.828 | 0.360 / 19.77 / 0.701 |
| vector mean→max | 0.323 / 14.61 / 0.829 | 0.303 / 15.97 / 0.809 | 0.380 / 18.76 / 0.703 |
| vector max→max | 0.335 / 14.20 / 0.872 | 0.285 / 11.47 / 0.866 | 0.381 / 20.56 / 0.762 |

**Interpretation.** Milestone 11’s gate: ATS scalar-vs-vector was within noise; nine-pack separates them. **Scalar max/max wins the two highest-weighted tasks.** Cold vector is not primary. See [`vector_vs_scalar.md`](vector_vs_scalar.md).

---

## 2. Age-covariates ablation *(rejected)*

| Metric | P2-G baseline | + tissue/sex conditioning | Δ |
|--------|-------------:|--------------------------:|--:|
| Tissue F1 | 0.355 | 0.340 | −0.015 |
| Age MAE | 13.431 | 14.135 | +0.704 |
| Sex AUROC | 0.853 | 0.834 | −0.019 |

**Interpretation.** Linear embedding-concat age head is the wrong lever. **Do not adopt.** Unconditioned P2-G stays the cascade finalist.

---

## 3. N-light — **`rho_hidden=64` default**

| Arm | folds | Tissue F1 | Age MAE | Sex AUROC |
|-----|------:|----------:|--------:|----------:|
| P2-G cascade | 3 | **0.355** | **13.431** | 0.853 |
| N-light rho=10 (historical) | 3 | 0.273 | 15.057 | 0.742 |
| **N-light rho=64 (default)** | 3 | **0.308** | **14.727** | **0.852** |

Canonical run ids: `stage0-7h-nine-pack-m-only-wide-f*`. Product configs + `flat_region` defaults updated 2026-09-08.

**Interpretation.** Wider rho **helps** one-hop (especially sex) but **does not close** the cascade gap (~5 F1 / ~1.3 y). Milestone **12** still pairs **P2-G + N-light@64**. Do not cite stale `external_test` 0.194/20.6 for this comparison.

---

## 4. Vector warm-start (LP-FT) — **max→max done; mean→max in flight**

Hypothesis: cold vector lost on optimization, not representation. Transplant scalar checkpoint → freeze ~4 ep → fine-tune at reduced LR.

### 4a. Vector max→max warm (complete, 3 folds)

| | Tissue F1 | Age MAE | Sex AUROC |
|--|----------:|--------:|----------:|
| Cold vector max→max | 0.333 | 15.408 | 0.834 |
| **Warm vector max→max** | **0.342** (±0.015) | **13.406** (±0.544) | **0.851** (±0.076) |
| P2-G scalar (ref) | 0.355 | 13.431 | 0.853 |

Folds (warm): f0 **0.348 / 12.88 / 0.949**, f1 0.321 / 13.18 / 0.840, f2 0.355 / 14.16 / 0.765.

**Interpretation.** Warm-start is a **real win vs cold vector** (esp. age + fold stability). Vs P2-G: age/sex essentially tied; tissue still **−0.013**. **Not enough to replace P2-G** as cascade finalist. Keep as diagnostic evidence that vector underperformance was partly optimization.

### 4b. Vector mean→max warm (in flight)

| | Tissue F1 | Age MAE | Sex AUROC |
|--|----------:|--------:|----------:|
| Cold vector mean→max | 0.335 | 16.448 | 0.780 |
| Warm fold 0 only | 0.364 | 12.799 | 0.941 |
| Warm folds 1–2 | **running** (~f1 ep12) | — | — |

**Interpretation.** Fold 0 matches P2-G fold 0 almost exactly — promising, but **one fold ≠ decision**. Finish 3-fold mean before any finalist change.

---

## 5. ATS seed-43 pooling (Track B.4) — done

| Arm | folds | Tissue F1 | Age MAE | Sex AUROC |
|-----|------:|----------:|--------:|----------:|
| P2-G-s2 | 3 | **0.379** | **19.251** | **0.765** |
| scalar-mean-max-s2 | 3 | 0.352 | 22.200 | 0.706 |
| scalar-max-mean-s2 | 3 | 0.353 | 20.824 | 0.703 |
| vector-mean-max-s2 | 3 | 0.342 | 21.735 | 0.658 |

**Interpretation.** Second seed still prefers P2-G. Do not mix ATS absolute age with nine-pack. Details: [`../stage0_7h_ats_pooling_s2/analysis.md`](../stage0_7h_ats_pooling_s2/analysis.md).

---

## 6. One-hop correctness smokes *(queued behind warm-start)*

Cheap checks (ATS fold 0, ~5 ep, **`rho_hidden=64`**):

1. Seed-gene mask **G0 vs G1** on `flat_region` (never tested on one-hop; 9c = cascade only).
2. Multi-seed restarts **42 / 43 / 44**.

**Status:** ImportError fixed (`open_betas_for_matrix` re-export + `sample_ids=`). Waiter
`scratch/logs/7h_onehop_correctness.*` starts when `run_7h_age_cov_and_warmstart_queue.sh`
exits. Report → [`../stage0_7h_onehop_correctness/`](../stage0_7h_onehop_correctness/).

**Interpretation.** Still the right **validate-before-scale** gate for N-light tricks /
Milestone 12 light OOF. Until G0/G1 + multi-seed diversity land, do not assume those
paths work on one-hop.

---

## 7. Trait expansion hygiene (10c)

See [`trait_hygiene.md`](trait_hygiene.md), [`label_status_census.md`](label_status_census.md).

| Trait | Status |
|-------|--------|
| age / sex / tissue | wired (tissue CE not expanded) |
| disease / cancer | `label_status` ready (~12.2k / ~9.1k case+control); GPU freeze-reuse next |
| bmi / ancestry | label-prep / heads stubbed in tree; GPU gated (≥1k bar / review) |
| blood / brain | **defer** heads |

**Interpretation.** Live train configs still **age/tissue/sex only**. No GPU trait expansion until one-hop smokes + freeze-reuse plan are green.

---

## 8. Milestone 12 skeleton numbers

| Finalist | Tissue F1 | Age MAE | Sex AUROC |
|----------|----------:|--------:|----------:|
| P2-G cascade | **0.355** | **13.431** | 0.853 |
| N-light@64 | 0.308 | 14.727 | 0.852 |

**Interpretation.** OOF = these two only. rho10 m-only (0.273) is historical, not the light deliverable.

---

## 9. Loader self-check

```json
{"matrix_id": "matrix-hub-nine-pack-virtual-v1", "shape": [34234, 482379], "checked_rows": 8, "checked_cols": 64, "finite_fraction": 0.91796875, "packs_in_check": {"matrix-hub-blood-full-v1": 8}, "platform_claim": "HM450_only_no_cross_platform"}
```

**Interpretation.** HM450-only; no cross-platform claim.

---

## 10. Next steps (ordered)

1. **Finish vector mean→max warm-start** (folds 1–2) → update §4b mean; if 3-fold mean still ≤ P2-G tissue, **close warm-start** without changing finalists.
2. **Run one-hop correctness smokes** (auto after warm queue, or launch manually on GPU 0 when free) → seed-mask + multi-seed report.
3. **Only then:** Milestone **12** scheduling (5×6 OOF on P2-G + N-light@64) and/or disease/cancer **freeze-and-reuse** probes (no joint retrain).
4. Soft-stop: do **not** auto-launch Stage B GPU / full OOF / pack-mask trait heads.
5. Optional follow-up only if warm vector still lags: stage-1 `region_pool: mean` pretrain-then-transplant (known dense-grad issue under max pooling).

**Does not authorize yet:** vector as primary; OOF on all arms; blood/brain heads; trusting one-hop seed-mask/multi-seed without the smoke report.
