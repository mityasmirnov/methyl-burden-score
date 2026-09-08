# Milestone 10 — nine-pack campaign analysis

Updated: `2026-09-08T15:10+02:00`

- Matrix: `matrix-hub-nine-pack-virtual-v1` · Screen split: `hub-nine-pack-3fold-v1` (**34 234**)
- OOF split: `hub-nine-pack-5fold-v1` · Platform: **HM450 this OOF**; encoder is gene-invariant for later EPIC/ONT
- **Product readout (locked):** frozen **`mbs_enet_nested` / `rbs_enet`**; `mbs_e2e` is encoder diagnostic only
- Detail dump: [`vector_vs_scalar.md`](vector_vs_scalar.md)
- Pre-OOF recipe: [`docs/plans/milestone-10e-staged-rbs-mbs-training.md`](../../../docs/plans/milestone-10e-staged-rbs-mbs-training.md)

---

## Live status

| Job | Status | Detail |
|-----|--------|--------|
| Vector warm-starts (max + mean) | **done** | See leaderboard |
| **One-hop smokes** | **done** | G0 ≫ G1; seeds 42/43/44 distinct (span 0.027) |
| Dense stage-1 (S1) | GPU 2 **fold 1 ~epoch 6/15** | then **hand off to N-light OOF** (skip transplants) |
| N-light nested enet (3-fold) | **done** | **0.368 / 9.88 / 0.803** vs e2e 0.308 / 14.73 / 0.852 |
| **12 N-light 5×6** | **queued GPU 2** | after S1; max VRAM |
| Cascade 10e S1–S4 | pending | still gates cascade 5×6 |
| Milestone 11 Stage B | deferred parallel | does **not** block OOF |

---

## Executive decisions

| Decision | Status |
|----------|--------|
| Cascade *topology* | **P2-G** (scalar max/max still the best *joint-e2e* arm) |
| How to *train* it before OOF | **Staged S1–S4** (dense vector RBS → freeze → MBS hop → unfreeze). Do **not** 5×6 joint e2e as-is |
| Light *encoder* | **N-light@64** |
| Light *product readout* | **`mbs_enet_nested` 3-fold 0.368 / 9.88 / 0.803** — no longer trails cascade e2e |
| Milestone **12** | **N-light 5×6 first**; cascade still gated on 10e; **not gated on 11** |
| Cold / warm vector | Warm closes gap; **neither replaces P2-G topology** |
| Age covariates | **Rejected** |

---

## Why joint `mbs_e2e` is the wrong OOF recipe

The campaign ranked arms with jointly trained multitask heads. That ranking is **not** the best predictor these encoders produce.

1. **Frozen RBS + classical beats gene MBS e2e.** Best campaign tissue F1 **0.368** and best age MAE **9.638** are `rbs_linear_probe`, not any e2e arm. P2-G: e2e 0.355 / 13.43 vs RBS probe **0.364 / 12.56**; MBS linear probe age **11.57**.
2. **LP-FT already showed freeze-then-train works.** Cold vector lost; warm-start from a frozen scalar encoder closed most of the gap. Jointly training a random gene hop with the encoder was the same class of mistake as joint e2e heads.
3. **enet ≫ e2e, especially N-light.** Nine-pack N-light@64 3-fold nested:
   **0.368 / 9.88 / 0.803** vs e2e **0.308 / 14.73 / 0.852**. Same MBS
   features: `mbs_linear_probe` also beats e2e — **readout, not encoder**.

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

Cheap **1-fold** S1–S4 smoke still gates **cascade** 5×6. **N-light 5×6 does not wait.** Dense stage-1 (GPU 2) is **S1 only**.

---

## Data: 34k encoder vs 173k freeze-reuse

Nine-pack **34 234** (HM450) is the honest **encoder-training** cohort. Catalog Hub **173 076** / 1 763 studies is **frozen-score reuse**, not a leaked joint retrain.

| Trait | n (honest labels) | How |
|-------|------------------|-----|
| ATS / age / tissue / sex | nine-pack split | encoder S1–S4 + frozen enet |
| BMI | 2 070 | freeze-reuse regression |
| Ancestry | 1 380; no class ≥1k | freeze-reuse CE; weak per-class |
| Cancer types | pack ~9.1k case+control; diagnosis strings | `label_status` + `phenotype_value`; **not** `cancer_mask` |
| Individual diseases | n_cases ≥ 600 | **Alzheimer’s 945 only**; schizophrenia 536 near-miss |
| Blood / brain packs | catalogue / control-only | **no trait head** |

Hub disease-pack cases: Alzheimer’s **945** (≥600), schizophrenia **536**, SLE **341**, Parkinson’s **333**, UC **258**, MS **228**, RA **225**, psoriasis **211**, stroke **204**, … (28 labels). Freeze-reuse only n_cases≥600 + matched `label_status=control`.

---

## Results digest (vector_vs_scalar)

1. **Cold 5-combo:** P2-G best joint-e2e tissue+age; cold vector loses.
2. **rho 10→64:** adopted for N-light.
3. **Age covariates:** rejected.
4. **Warm-starts:** optimization helps; no topology flip.
5. **RBS/MBS probes beat e2e** — now the **training/readout lock**, not a side note.

---

## One-hop smokes — **done**

G0 (dense): **0.335 / 17.3 / 0.750**. G1 seed-mask **collapses** (0.146 / 21.7 / 0.56) — **not adopted**. Multi-seed 42/43/44 distinct (tissue F1 span 0.027).  
Report: `reports/inspection/stage0_7h_onehop_correctness/`.

---

## Outlook / next steps

1. GPU 2: finish dense S1, **skip transplants**, start **N-light 5×6**.
2. CPU: finish P2-G nested enet 3-fold; re-rank under enet.
3. Freeze-reuse n≥600 (AD 945, BMI, ancestry, cancer pack).
4. Cascade 5×6 only after 10e S1–S4 smoke.
5. Keep encoder gene-invariant for later EPIC/ONT; this OOF stays HM450.

---

## Trait hygiene (10c)

age/sex/tissue wired · disease/cancer `label_status` ready · BMI/ancestry joined/stubbed · blood/brain defer.  
See [`trait_hygiene.md`](trait_hygiene.md).
