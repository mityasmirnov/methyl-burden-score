# Milestone 10 — nine-pack campaign analysis

Updated: `2026-09-08` (manual synthesis over auto `summary.json` / `vector_vs_scalar.*`)

- Matrix: `matrix-hub-nine-pack-virtual-v1` (virtual multi-store)
- Split: `hub-nine-pack-3fold-v1` (**34 234** samples)
- Platform: **HM450 only** — no cross-platform claim
- Primary readout: **`mbs_e2e`** on outer **test** (means ± population SD across 3 folds unless noted)

**How to read this file:** newest / decision-blocking results first; each section ends with an explicit **Interpretation**. Machine snapshots that only cover P2-G + m-only live in [`summary.json`](summary.json); the pooling grid dump is [`vector_vs_scalar.md`](vector_vs_scalar.md).

---

## Executive decisions (current)

| Decision | Status | Basis |
|----------|--------|-------|
| Cascade finalist | **P2-G scalar max/max** | Best tissue F1 **and** best age MAE in the full 5-combo nine-pack pooling grid |
| Light-model finalist | **N-light m-only (`rho_hidden=64`)** | Wide screen beat rho=10 on tissue/age/sex; adopted as default 2026-09-08 |
| Vector RBS at nine-pack | **Does not beat scalar P2-G** | Both vector arms trail on tissue + age |
| Age head + tissue/sex covariates | **Not adopted** | All three headline metrics slightly worse |
| Widen m-only `rho_hidden` 10→64 | **Adopted as N-light default** | 0.273→**0.308** tissue, 15.1→**14.7** age, 0.742→**0.852** sex; still trails cascade |
| Milestone **12** OOF scope | **Finalists only** | P2-G cascade + N-light — not all pooling arms, not all 9 trait packs |
| Disease / cancer / blood / brain heads | **Blocked / deferred** | Need `label_status` wiring; blood/brain pack semantics broken |

**Still in flight / blocked**

- Vector max/max **warm-start from P2-G** — fold 0 done; fold 1 training (do not treat as final).
- One-hop correctness smokes (seed-mask + multi-seed) — **blocked** on `ImportError: open_betas_for_matrix`.
- Trait GPU expansion — not started; age/tissue/sex heads only in live configs.

---

## 1. Scalar vs vector RBS at nine-pack scale *(latest scale screen)*

Matched budget: **15 epochs × 3 folds**, same split/seed/checkpoint rule as P2-G.

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

**Interpretation.** This is the gate Milestone 11 was waiting on: ATS-scale scalar-vs-vector was within noise; nine-pack (~34k) is large enough to separate them. **Scalar max/max wins on the two highest-weighted tasks** (tissue weight 3.0, age primary). Vector embeddings do not buy a tissue or age win here — both vector arms sit ~2 F1 points and ~2–3 years of MAE behind P2-G. Scalar mean/max is the only arm that beats P2-G on sex (+0.022 AUROC) while staying close on age, but it loses tissue; with current loss weights that is not enough to dethrone P2-G. **Cascade finalist is no longer provisional: lock P2-G scalar max/max** for Stage B / OOF planning. See also [`vector_vs_scalar.md`](vector_vs_scalar.md).

---

## 2. Age-covariates ablation *(tissue + sex → age head)*

Full 3-fold rerun of P2-G with `model.age_covariates: [tissue, sex]` (ground-truth embeddings concatenated into the age head only).

| Metric | P2-G baseline | + tissue/sex conditioning | Δ |
|--------|-------------:|--------------------------:|--:|
| Tissue F1 | 0.355 | 0.340 | −0.015 |
| Age MAE | 13.431 | 14.135 | +0.704 (worse) |
| Sex AUROC | 0.853 | 0.834 | −0.019 |

Per-fold age MAE: baseline `[12.80, 16.00, 11.49]` vs conditioned `[12.75, 13.27, 16.38]` — fold 2 flips from best→worst and dominates the mean.

**Interpretation.** The biological motivation (clocks are tissue-dependent; blood easy / brain+colon hard) is real, but a linear embedding-concat age head is the wrong lever on this encoder: headlines move **against** conditioning, and the blood-vs-brain median gap barely narrows. **Do not adopt.** Do not spend on FiLM/nonlinear variants without new evidence. Unconditioned P2-G remains the cascade finalist.

---

## 3. N-light (one-hop m-only) — **`rho_hidden=64` is now the default**

| Arm | folds | Tissue F1 | Age MAE | Sex AUROC |
|-----|------:|----------:|--------:|----------:|
| P2-G cascade | 3 | **0.355** | **13.431** | 0.853 |
| N-light m-only (`rho_hidden=10`, historical) | 3 | 0.273 (±0.050) | 15.057 (±2.150) | 0.742 (±0.110) |
| **N-light m-only (`rho_hidden=64`, default)** | 3 | **0.308** (±0.053) | **14.727** (±1.370) | **0.852** (±0.032) |

Configs updated 2026-09-08: `stage0_7h_nine_pack_m_only.yaml`, light-mean/max product
YAMLs, and `flat_region` code defaults. Canonical nine-pack numbers for OOF planning
are the **rho=64** row (`stage0-7h-nine-pack-m-only-wide-f*`).

**Interpretation.** At ATS scale, N-light was near-parity with cascade; at nine-pack
scale a gap remains even at rho=64 (~5 tissue F1 points, ~1.3 y age). Widening the
DeepSet rho MLP from 10→64 **clearly helps** vs the narrow bottleneck (especially
sex AUROC, which catches cascade), so **64 is the N-light default going forward**.
It still does **not** replace cascade — Milestone **12** keeps both finalists
(P2-G + N-light@64).

---

## 4. Vector warm-start from P2-G *(in flight)*

Hypothesis: cold vector underperforms because of harder optimization, not because vector RBS is inherently worse — warm-start encoder/heads from scalar P2-G.

| Status | Tissue F1 | Age MAE | Sex AUROC |
|--------|----------:|--------:|----------:|
| Cold vector max→max (3-fold mean) | 0.333 | 15.408 | 0.834 |
| Warm-start fold 0 only | 0.348 | 12.884 | 0.949 |
| Warm-start folds 1–2 | **running** | — | — |

**Interpretation.** Fold 0 looks encouraging vs cold vector on the same fold (0.335 / 14.20 / 0.872 → 0.348 / 12.88 / 0.949), but **one fold is not a decision**. Wait for the full 3-fold mean before any claim that warm-start resurrects vector RBS — and even a warm-start win would need to beat **scalar P2-G** (0.355 / 13.43), not merely beat cold vector.

---

## 5. ATS seed-43 pooling (Track B.4) — secondary confirmation

Seed **43** on `hub-ats-7e-3fold-v1` (exercises the seed-offset fix). Full write-up: [`../stage0_7h_ats_pooling_s2/analysis.md`](../stage0_7h_ats_pooling_s2/analysis.md).

| Arm | folds | Tissue F1 | Age MAE | Sex AUROC |
|-----|------:|----------:|--------:|----------:|
| P2-G-s2 | 3 | **0.379** (±0.054) | **19.251** | **0.765** |
| scalar-mean-max-s2 | 3 | 0.352 | 22.200 | 0.706 |
| scalar-max-mean-s2 | 3 | 0.353 | 20.824 | 0.703 |
| vector-mean-max-s2 | 3 | 0.342 | 21.735 | 0.658 |

**Interpretation.** On a second seed, **P2-G still leads** every column in this incomplete 2×2 (+P2) grid. Useful as a diversity / correctness check after the seed-offset bugfix; **not** a reason to reopen the nine-pack pooling decision. ATS age MAEs are on a different scale/cohort than nine-pack — do not mix absolute age numbers across these tables.

---

## 6. One-hop correctness smokes *(fix applied — rerunning)*

Intended cheap checks (fold 0, ~5 ep) on `flat_region` with **`rho_hidden=64`**:

1. Seed-gene mask G0 vs G1 (9c only ever tested this on **cascade**).
2. Multi-seed restarts 42/43/44 (path unproven on one-hop after the seed-offset fix).

**Status:** queue previously hit `ImportError: open_betas_for_matrix` (wrong import
from `mbs.matrix.store`). **Fixed:** re-export added on `store.py`, smoke script
uses keyword `sample_ids=…`, and training uses rho=64. Report lands in
[`../stage0_7h_onehop_correctness/`](../stage0_7h_onehop_correctness/).

**Interpretation.** These remain the right **validate-before-scale** steps before
any multi-trait / multi-restart N-light campaign. Do **not** launch Milestone 12
light-model OOF prep that assumes seed-mask or multi-seed work on one-hop until
this report shows stable G0/G1 runs and non-zero multi-seed diversity.

---

## 7. Trait expansion hygiene (10c) — data, not modeling

Policy + census: [`trait_hygiene.md`](trait_hygiene.md), [`label_status_census.md`](label_status_census.md). Sidecar parquet written; **live training heads remain age / tissue / sex only**.

| Trait family | Trainable now? | Note |
|--------------|----------------|------|
| age / sex | yes | Wired in P2-G / m-only configs |
| tissue | yes, constrained | Only `whole blood` clears ≥1k among 64 labels — do **not** expand CE head |
| disease | labels ready | `label_status` case+control ≈ **12.2k**; must not use pack `disease_mask` as control |
| cancer | labels ready | case+control ≈ **9.1k**; same mask caveat |
| blood / brain | **defer** | Pack rows are catalogue `sample_type=control`, not case/control |
| bmi / ancestry | **not in table** | Separate small files; not joined |

**Interpretation.** “All 9 packs” is not an executable training plan today. The honest multi-trait light-model target is **≤5** traits after wiring disease/cancer heads to `label_status`, and only after one-hop smokes pass. Training on naive pack masks would silently reproduce the case≠control bug.

---

## 8. Baseline reference pair (P2-G vs m-only) — auto snapshot

Generated machine view (may lag the sections above if only these two arms were refreshed): [`summary.json`](summary.json).

| Arm | folds | Tissue F1 | Age MAE | Sex AUROC |
|-----|------:|----------:|--------:|----------:|
| P2-G cascade | 3 | 0.355 (±0.054) | 13.431 (±2.319) | 0.853 (±0.087) |
| one-hop m-only | 3 | 0.273 (±0.061) | 15.057 (±2.634) | 0.742 (±0.135) |

### Fold detail

```json
{
  "P2-G": [
    {"fold": 0, "tissue_f1": 0.364, "age_mae": 12.796, "sex_auroc": 0.950, "best_epoch": 12},
    {"fold": 1, "tissue_f1": 0.297, "age_mae": 16.001, "sex_auroc": 0.830, "best_epoch": 4},
    {"fold": 2, "tissue_f1": 0.404, "age_mae": 11.494, "sex_auroc": 0.780, "best_epoch": 10}
  ],
  "m_only": [
    {"fold": 0, "tissue_f1": 0.343, "age_mae": 12.073, "sex_auroc": 0.890, "best_epoch": 13},
    {"fold": 1, "tissue_f1": 0.245, "age_mae": 17.055, "sex_auroc": 0.708, "best_epoch": 6},
    {"fold": 2, "tissue_f1": 0.230, "age_mae": 16.045, "sex_auroc": 0.626, "best_epoch": 7}
  ]
}
```

**Interpretation.** This pair is the **Milestone 12 OOF skeleton**: cascade product path + cheap light model. m-only fold 0 is competitive with P2-G on that fold alone; folds 1–2 drag the mean — another reason cheap multi-seed / seed-mask smokes on one-hop matter before spending OOF budget.

---

## 9. Loader self-check

```json
{"matrix_id": "matrix-hub-nine-pack-virtual-v1", "shape": [34234, 482379], "checked_rows": 8, "checked_cols": 64, "finite_fraction": 0.91796875, "packs_in_check": {"matrix-hub-blood-full-v1": 8}, "platform_claim": "HM450_only_no_cross_platform"}
```

**Interpretation.** Virtual multi-store loads; finite fraction ~0.92 on the probe slice is expected for sparse/missing methylation. **HM450-only** — do not generalize scores to EPIC/other arrays from these runs.

---

## 10. What this does *not* authorize

- Declaring vector RBS primary (warm-start unfinished; cold vector lost).
- Expanding disease/cancer GPU heads on pack-membership masks.
- Blood/brain/bmi/ancestry as “nine arms.”
- Milestone **12** OOF on the full Stage A matrix or all pooling combos.
- Trusting one-hop seed-mask / multi-seed behavior before the import fix + smokes.

**Next engineering priorities:** (1) finish vector warm-start 3-fold and record vs P2-G; (2) fix + run one-hop correctness smokes; (3) only then consider light-model multi-trait or Milestone 12 scheduling.
