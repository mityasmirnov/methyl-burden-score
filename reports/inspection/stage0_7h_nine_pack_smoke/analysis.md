# Milestone 10 — nine-pack campaign analysis

Updated: `2026-09-08T17:05+02:00`

- Matrix: `matrix-hub-nine-pack-virtual-v1` · Screen split: `hub-nine-pack-3fold-v1` (**34 234**)
- OOF split: `hub-nine-pack-5fold-v1` · Platform: **HM450 this OOF**; encoder is gene-invariant for later EPIC/ONT
- **Product readout (locked):** frozen **`mbs_enet_nested` / `rbs_enet`**; `mbs_e2e` is encoder diagnostic only
- Encoder aux: disease classes + cancer types with nine-pack disease-tissue **n>200**
- Detail dump: [`vector_vs_scalar.md`](vector_vs_scalar.md)
- Pre-OOF recipe: [`docs/plans/milestone-10e-staged-rbs-mbs-training.md`](../../../docs/plans/milestone-10e-staged-rbs-mbs-training.md)

---

## Live status

| Job | Status | Detail |
|-----|--------|--------|
| Vector warm-starts | **done** | See leaderboard |
| One-hop smokes | **done** | G0 ≫ G1; seeds 42/43/44 distinct |
| Dense S1 (mean/mean) | **fold 0 done**; fold 1 stopped ~epoch 12/15 | GPU 2 handed to N-light OOF. Fold 0 e2e **0.348 / 14.52**; linear age **10.93** |
| N-light nested enet (3-fold screen) | **done** | **0.368 / 9.88 / 0.803** vs e2e 0.308 / 14.73 / 0.852 |
| P2-G nested enet | **3 / 3 folds** | MBS **0.335 / 9.81 / 0.759**; RBS age MAE **19.8** (overfits) |
| Freeze-reuse P2-G | **done** | Pack cancer **0.954**; Alzheimer’s **0.838**; cancer subtype macro-F1 **0.202**; BMI not useful |
| **12 N-light 5×6** | **running GPU 2** | **30 ep** + n>200 aux; `f0-r0` **ep ~19–20/30** (1/30). Val disease ~0.76 / cancer ~0.85. 16-ep archived **0.300 / 8.63 / 0.880** |
| 10e S1–S4 1-fold smoke | **done (negative)** | fold-0 e2e **0.300 / 14.97 / 0.900** loses to native P2-G **0.364 / 12.80 / 0.950** — no 3-fold staged recipe |
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
| Disease classes n>200 | 9 diagnoses (AD 945 … childhood asthma 202) | **encoder aux BCE** (Milestone 12) |
| Cancer types n>200 | 9 types (glioma 481 … sarcoma 239) | **encoder aux BCE** (Milestone 12) |
| BMI | 2 070 | freeze-reuse regression |
| Ancestry | 1 380; no class ≥1k | freeze-reuse CE; weak per-class |
| Remaining diagnoses / types | n≤200 or exact 200 | freeze-reuse / census; **not** pack masks |
| Blood / brain packs | catalogue / control-only | **no trait head** |

Hub disease-pack cases in nine-pack *disease tissue*: Alzheimer’s **945**, schizophrenia **536**, SLE **341**, Parkinson’s **333**, RA **250**, MS **228**, psoriasis **211**, stroke **204**, childhood asthma **202** (encoder aux). Exact-200 and below stay out of joint training. Pack-matched `control` = negatives; adjacent-normal = unknown.

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

1. Let **N-light 5×6** finish on GPU 2 (30-ep + n>200 aux; `f0-r0` ~ep 20/30).
2. P2-G nested **3/3 done** (MBS **0.335 / 9.81 / 0.759** vs N-light 3-fold **0.368 / 9.88 / 0.803**).
3. **10e staged S1–S4 smoke negative** — do **not** regenerate S1 folds 1–2. Cascade 5×6, if launched, uses **native P2-G**, not the staged recipe.
4. Alzheimer’s freeze-reuse **0.838** (3-fold); cancer subtypes distinguishable (macro-F1 0.202). Pack cancer still strong.
5. Keep encoder gene-invariant for later EPIC/ONT; this OOF stays HM450.

---

## Trait hygiene (10c)

age/sex/tissue wired · disease/cancer `label_status` ready · BMI/ancestry joined/stubbed · blood/brain defer.  
See [`trait_hygiene.md`](trait_hygiene.md) and
[`disease_cancer_frozen_probes/analysis.md`](disease_cancer_frozen_probes/analysis.md)
(cancer AUROC ~0.95; broad disease ~0.59; BMI not useful).
