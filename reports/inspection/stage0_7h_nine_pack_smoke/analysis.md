# Milestone 10 — nine-pack campaign analysis

Updated: `2026-09-08T14:40+02:00`

- Matrix: `matrix-hub-nine-pack-virtual-v1` · Split: `hub-nine-pack-3fold-v1` (**34 234**)
- Platform: **HM450 only** · Architecture grid primary: **`mbs_e2e`** outer test
- **Product readout (locked for planning):** frozen **`rbs_enet` / `mbs_enet(_nested)`** co-primary; `mbs_e2e` secondary
- Detail dump: [`vector_vs_scalar.md`](vector_vs_scalar.md)
- Pre-OOF recipe: [`docs/plans/milestone-10e-staged-rbs-mbs-training.md`](../../../docs/plans/milestone-10e-staged-rbs-mbs-training.md)

---

## Live status

| Job | Status | Detail |
|-----|--------|--------|
| Vector warm-starts (max + mean) | **done** | See leaderboard |
| **One-hop smokes** | **running** GPU 0 | G0 ≈0.335 / 17.3 / 0.75; G1 in flight |
| Dense stage-1 (S1 experiment) | other GPU | informs 10e; **not** a 5×6 launch |
| **10e staged recipe smoke** | **pending — blocks 12** | S1–S4 + enet; not joint e2e |
| Milestone 11 Stage B | deferred parallel | does **not** block OOF |

---

## Executive decisions

| Decision | Status |
|----------|--------|
| Cascade *topology* | **P2-G** (scalar max/max still the best *joint-e2e* arm) |
| How to *train* it before OOF | **Staged S1–S4** (dense vector RBS → freeze → MBS hop → unfreeze). Do **not** 5×6 joint e2e as-is |
| Light *encoder* | **N-light@64** |
| Light *product readout* | **`mbs_enet_nested` required** (ATS: e2e age 17 vs nested enet 10) |
| Cold / warm vector | Warm closes gap; **neither replaces P2-G topology** |
| Age covariates | **Rejected** |
| Milestone **12** | Blocked on **10e** recipe smoke + enet co-primary; **not gated on 11** |

---

## Why joint `mbs_e2e` is the wrong OOF recipe

The campaign ranked arms with jointly trained multitask heads. That ranking is **not** the best predictor these encoders produce.

1. **Frozen RBS + classical beats gene MBS e2e.** Best campaign tissue F1 **0.368** and best age MAE **9.638** are `rbs_linear_probe`, not any e2e arm. P2-G: e2e 0.355 / 13.43 vs RBS probe **0.364 / 12.56**; MBS linear probe age **11.57**.
2. **LP-FT already showed freeze-then-train works.** Cold vector lost; warm-start from a frozen scalar encoder closed most of the gap. Jointly training a random gene hop with the encoder was the same class of mistake as joint e2e heads.
3. **enet ≫ e2e, especially N-light.** ATS N-light: `mbs_e2e` ~0.378 / age **17.1** vs `mbs_enet_nested` **0.387 / 10.3**. ATS P2-G: `mbs_enet` 0.385 vs `mbs_e2e` 0.373. Nine-pack N-light@64 still has **no** enet score (`include_mbs_enet` was off) — CPU before GPU OOF.

**Interpretation.** Region RBS carries more signal than max-pooled gene MBS for several tasks; jointly trained heads throw some of it away. Train the encoder in stages; **report frozen elastic-net**. Improving e2e (deeper / separate heads) is a head problem, not a reason to skip enet or to launch 5×6 on today’s joint recipe.

---

## Leaderboard (`mbs_e2e`) — architecture screen only

| Arm | Tissue F1 | Age MAE | Sex AUROC |
|-----|----------:|--------:|----------:|
| **P2-G scalar max/max** | **0.355** | **13.431** | 0.853 |
| vector mean→max **warm** | 0.345 | 14.578 | 0.862 |
| vector max→max **warm** | 0.342 | **13.406** | 0.851 |
| scalar mean/max | 0.330 | 13.496 | **0.875** |
| vector mean→max cold | 0.335 | 16.448 | 0.780 |
| vector max→max cold | 0.333 | 15.408 | 0.834 |
| scalar max/mean | 0.318 | 17.426 | 0.812 |
| **N-light@64** | 0.308 | 14.727 | 0.852 |
| N-light@10 | 0.273 | 15.057 | 0.742 |
| P2-G + age covariates | 0.340 | 14.135 | 0.834 |

This table chose the **topology**. It does **not** choose the OOF training recipe or the product readout.

---

## Staged training (10e) — lock before 12

| Stage | Train | Freeze | Pool |
|------:|-------|--------|------|
| **S1** | CpG→region (`phi`, `region_rho`), prefer **vector** RBS | — | **mean** (dense grads; max starves `region_rho`) |
| **S2** | Scalar / burden RBS (max or learned) | S1 frozen, then low-LR | product aggregation |
| **S3** | Gene MBS (`gene_rho` / learned pooler — **not max-only** if a small encoder wins the smoke) | RBS frozen | gene hop |
| **S4** | Joint fine-tune @ ~3e-4 | — | keep S2/S3 pooling |

Then freeze the encoder and fit **enet / logistic** on RBS and MBS for ATS and extra traits.

Cheap **1-fold** S1–S4 smoke must beat P2-G joint e2e **and** match that run’s own RBS/MBS enet before any 5×6. Dense stage-1 (GPU 2) is **S1 only**.

---

## Data: 34k encoder vs 173k freeze-reuse

Nine-pack **34 234** (HM450) is the honest **encoder-training** cohort. Catalog Hub **173 076** / 1 763 studies is **frozen-score reuse**, not a leaked joint retrain.

| Trait | n (honest labels) | How |
|-------|------------------|-----|
| ATS / age / tissue / sex | nine-pack split | encoder S1–S4 + frozen enet |
| BMI | 2 070 | freeze-reuse regression |
| Ancestry | 1 380; no class ≥1k | freeze-reuse CE; weak per-class |
| Cancer types | pack ~9.1k case+control; diagnosis strings | `label_status` + `phenotype_value`; **not** `cancer_mask` |
| Individual diseases | see census below | freeze-reuse logistic; skip if n≪1k |
| Blood / brain packs | catalogue / control-only | **no trait head** |

Hub disease-pack cases (nine-pack Hub, `sample_overview`): Alzheimer’s **945**, schizophrenia **536**, SLE **341**, Parkinson’s **333**, UC **258**, MS **228**, RA **225**, psoriasis **211**, stroke **204**, … (28 labels). **Census n + matched controls before GPU**; do not train on naive pack masks.

---

## Results digest (vector_vs_scalar)

1. **Cold 5-combo:** P2-G best joint-e2e tissue+age; cold vector loses.
2. **rho 10→64:** adopted for N-light.
3. **Age covariates:** rejected.
4. **Warm-starts:** optimization helps; no topology flip.
5. **RBS/MBS probes beat e2e** — now the **training/readout lock**, not a side note.

---

## One-hop smokes (in flight)

G0 (dense, no mask): tissue **0.335**, age **17.3**, sex **0.750** (5 ep ATS fold 0, rho=64).  
G1 (seed-mask) + multi-seed 42/43/44 next.  
Report: `reports/inspection/stage0_7h_onehop_correctness/`.

Correctness gate on the **encoder**, not a substitute for nested enet.

---

## Outlook / next steps

1. Finish one-hop G1 + multi-seed → publish smoke report.
2. **CPU:** enet on existing P2-G + N-light@64 checkpoints (nine-pack skipped `include_mbs_enet`).
3. **10e:** 1-fold S1–S4 smoke vs P2-G joint e2e and vs that smoke’s RBS/MBS enet.
4. Freeze-reuse: ATS, BMI, ancestry, cancer types, AD/PD/stroke/… on frozen scores.
5. **Then** Milestone **12** 5×6 on staged recipe + enet co-primary.
6. Milestone **11** stays parallel. Soft-stop: no joint-e2e 5×6, no pack-mask trait GPU, no auto Stage B.

---

## Trait hygiene (10c)

age/sex/tissue wired · disease/cancer `label_status` ready · BMI/ancestry joined/stubbed · blood/brain defer.  
See [`trait_hygiene.md`](trait_hygiene.md).
