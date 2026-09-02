# Stage 0 - 7G cascade tissue probe (P0-P5)

Frozen split `hub-ats-7e-3fold-v1`; **65536 loci / 1 restart**. P0-P4 use 15 epochs; P5 uses a 30-epoch ceiling with validation-tissue early stopping.

## Arm glossary

### Naming prefixes

| Prefix | Meaning | Examples |
|--------|---------|----------|
| **P** | **Probe** ablation (7G tissue investigation or 7G′ `-G` gene-only) | `P2-end2end-tissue-weight`, `P4-G` |
| **C-mvalue-** | **Classical** sklearn/HGB on M-values (no neural encoder) | `C-mvalue-enet`, `C-mvalue-enet-G` |
| **N-** | **Neural** encoder (CascadeDeepSet or flat/hier variants) | `N-cascade-l1`, `N-light-type` |
| **T-** | **Transparent** fixed-feature linear baselines | `T-mean-region`, `T-mean-gene` |
| **-G** | **Gene-linked** CpG panel only (`gene_cols`; 7G′ Stage A) | `P2-G`, `C-mvalue-ridge-G` |
| **-S** | Fold-safe **Selected** sparse panel (7G′ Stage B) | `C-mvalue-enetS`, `N-cascade-S` |

| Arm | Description |
|-----|-------------|
| `P0-baseline` | Replay 7G product cascade (`N-cascade-l1`): train + **late fusion** `[orphan_rbs | mbs | direct_contrib]` on all 65 536 prefix loci. |
| `P1-fusion-tissue-heavy` | **Fusion-only** on saved P0 scores: balanced logistic tissue head + PCA(32); no retrain. |
| `P2-end2end-tissue-weight` | Cascade **retrain** with tissue_loss_weight=3, age_loss_weight=0.3; test via **late fusion** (not MBS-only). |
| `P2-fusion-balanced` | P2 saved scores + balanced logistic late fusion (no retrain). |
| `P2-fusion-logistic` | P2 saved scores + multinomial logistic late fusion (no retrain). |
| `P2-fusion-sgd` | P2 saved scores + SGD logistic late fusion (no retrain). |
| `P3-region-head-bypass` | Transparent **`T-mean-region`**: region-mean M-values → linear heads; no CascadeDeepSet training. |
| `P4-fusion-balanced` | P4 saved scores + balanced logistic late fusion. |
| `P4-fusion-logistic` | P4 saved scores + logistic late fusion. |
| `P4-fusion-sgd` | P4 saved scores + SGD late fusion. |
| `P4-pooling-mean` | Like P2 loss weights but **mean/mean** CpG→region and region→gene pooling; late-fusion test. |
| `P5-epochs-30` | P2 loss weights, **max/max** pooling, **30-epoch** ceiling + early stop on validation tissue F1; late-fusion test. |
| `P5-mean-epochs-30` | P5 schedule with **mean/mean** pooling. |
| `C-mvalue-enet` | All prefix loci → M-value → elastic-net regression (age) / logistic enet (tissue, sex); **7G bake-off tissue winner** (F1≈0.334 on 65k). |
| `C-mvalue-enet-G` | Elastic-net on **identical gene-linked CpGs** as neural `-G` arms (fair Stage A comparator). |
| `C-mvalue-enetS` | **Stage B:** fold-safe stability-selected sparse CpG panel (outer-train only) → refit enet; one test eval per fold. |
| `P2-G` | CascadeDeepSet on **gene-linked CpGs only** (`gene_cols`); max/max pooling; 15 epochs; P2 loss weights; primary metric **`mbs_e2e`** (end-to-end MBS heads). |
| `P4-G` | Like P2-G but **mean/mean** pooling; same gene-only panel and **`mbs_e2e`** metric. |
| `N-cascade-l1` | **7F/7G product topology:** CascadeDeepSet (CpG→RBS→MBS) + fold-fitted direct enet + **late fusion** on saved blocks; Level-1 z on direct branch. |

## Comparable ranking (panel × eval mode)

Compare **only within the same row group** (same panel and eval mode). Stage A primary metric is **`mbs_e2e`** on the **gene-linked** panel; 7G tissue probe P0–P5 used **late fusion (`fusion_full`)** on the **65k prefix**. The same cascade checkpoint can score very differently under `mbs_e2e` vs `fusion_full`.

| Arm | Panel | Eval mode | Tissue macro-F1 | folds |
|-----|-------|-----------|----------------:|------:|
| `P0-baseline` | 65k prefix | `fusion_full` | 0.093 (±0.014) | 3 |
| `P0-baseline` | 65k prefix | `mbs_e2e` | — | 0 |
| `P2-end2end-tissue-weight` | 65k prefix | `fusion_full` | 0.376 (±0.059) | 3 |
| `P2-end2end-tissue-weight` | 65k prefix | `mbs_e2e` | — | 0 |
| `P4-pooling-mean` | 65k prefix | `fusion_full` | 0.360 (±0.053) | 3 |
| `P4-pooling-mean` | 65k prefix | `mbs_e2e` | — | 0 |
| `P5-epochs-30` | 65k prefix | `fusion_full` | 0.356 (±0.050) | 3 |
| `P5-epochs-30` | 65k prefix | `mbs_e2e` | 0.674 (±0.040) | 3 |
| `P2-G` | gene-linked | `fusion_full` | 0.380 (±0.051) | 3 |
| `P2-G` | gene-linked | `mbs_e2e` | 0.682 (±0.030) | 3 |
| `P4-G` | gene-linked | `fusion_full` | 0.374 (±0.048) | 3 |
| `P4-G` | gene-linked | `mbs_e2e` | 0.695 (±0.054) | 3 |
| `P5-G-mean` | gene-linked | `mbs_e2e` | 0.672 (±0.029) | 3 |
| `C-mvalue-enet` | 65k prefix | `classical` | 0.334 (±0.031) | 3 |
| `C-mvalue-enet-G` | gene-linked | `classical` | — | 0 |

**Fair pairs (examples):**

- **Late fusion, 65k:** `P2-end2end` vs `P4-pooling-mean` vs `C-mvalue-enet`.
- **Late fusion, gene-linked:** `P2-G` `fusion_full` ≈ `P2-end2end` on 65k (same loss weights).
- **MBS e2e, gene-linked:** `P2-G` vs `P4-G` vs `C-mvalue-enet-G` (Stage A lock metric).
- **MBS e2e, 65k:** backfill `evaluations.mbs_e2e` on P2/P4/P0 (P5 already has both).

## Per-arm summary

| Arm | Tissue macro-F1 | Balanced acc | Sex AUROC | Age MAE | Age R² |
|-----|-----------------|--------------|-----------|---------|--------|
| P0-baseline | 0.093 | 0.125 | 0.911 | 8.441 | 0.782 |
| P1-fusion-tissue-heavy | 0.124 | 0.131 | 0.911 | 8.441 | 0.782 |
| P2-end2end-tissue-weight | 0.376 | 0.422 | 0.885 | 8.682 | 0.767 |
| P2-fusion-balanced | 0.385 | 0.428 | 0.885 | 8.682 | 0.767 |
| P2-fusion-logistic | 0.376 | 0.422 | 0.885 | 8.682 | 0.767 |
| P2-fusion-sgd | 0.328 | 0.366 | 0.885 | 8.682 | 0.767 |
| P3-region-head-bypass | 0.389 | 0.424 | 0.908 | 7.946 | 0.802 |
| P4-fusion-balanced | 0.375 | 0.414 | 0.899 | 8.620 | 0.772 |
| P4-fusion-logistic | 0.360 | 0.404 | 0.899 | 8.620 | 0.772 |
| P4-fusion-sgd | 0.308 | 0.356 | 0.899 | 8.620 | 0.772 |
| P4-pooling-mean | 0.360 | 0.404 | 0.899 | 8.620 | 0.772 |
| P5-epochs-30 | 0.356 | 0.404 | 0.883 | 8.741 | 0.764 |
| P5-mean-epochs-30 | 0.361 | 0.405 | 0.893 | 8.597 | 0.775 |

## Locked 7G tissue comparator

`C-mvalue-enet`: mean F1 0.334. It remains the 7G bake-off reference on the 65k prefix until **C-mvalue-enet-G** (gene-linked panel) and **7G′ Stage B** (`C-mvalue-enetS`) complete.

## Per-fold tissue macro-F1

### P0-baseline

_Replay 7G product cascade (`N-cascade-l1`): train + **late fusion** `[orphan_rbs | mbs | direct_contrib]` on all 65 536 prefix loci._

- fold 0: F1=0.076, sex AUROC=0.806, age MAE=9.653
- fold 1: F1=0.104, sex AUROC=0.982, age MAE=7.566
- fold 2: F1=0.098, sex AUROC=0.943, age MAE=8.105

### P1-fusion-tissue-heavy

_**Fusion-only** on saved P0 scores: balanced logistic tissue head + PCA(32); no retrain._

- fold 0: F1=0.119, sex AUROC=0.806, age MAE=9.653
- fold 1: F1=0.119, sex AUROC=0.982, age MAE=7.566
- fold 2: F1=0.134, sex AUROC=0.943, age MAE=8.105

### P2-end2end-tissue-weight

_Cascade **retrain** with tissue_loss_weight=3, age_loss_weight=0.3; test via **late fusion** (not MBS-only)._

- fold 0: F1=0.308, sex AUROC=0.818, age MAE=9.923
- fold 1: F1=0.407, sex AUROC=0.932, age MAE=7.828
- fold 2: F1=0.412, sex AUROC=0.907, age MAE=8.295

### P2-fusion-balanced

_P2 saved scores + balanced logistic late fusion (no retrain)._

- fold 0: F1=0.309, sex AUROC=0.818, age MAE=9.923
- fold 1: F1=0.409, sex AUROC=0.932, age MAE=7.828
- fold 2: F1=0.437, sex AUROC=0.907, age MAE=8.295

### P2-fusion-logistic

_P2 saved scores + multinomial logistic late fusion (no retrain)._

- fold 0: F1=0.308, sex AUROC=0.818, age MAE=9.923
- fold 1: F1=0.407, sex AUROC=0.932, age MAE=7.828
- fold 2: F1=0.412, sex AUROC=0.907, age MAE=8.295

### P2-fusion-sgd

_P2 saved scores + SGD logistic late fusion (no retrain)._

- fold 0: F1=0.294, sex AUROC=0.818, age MAE=9.923
- fold 1: F1=0.332, sex AUROC=0.932, age MAE=7.828
- fold 2: F1=0.358, sex AUROC=0.907, age MAE=8.295

### P3-region-head-bypass

_Transparent **`T-mean-region`**: region-mean M-values → linear heads; no CascadeDeepSet training._

- fold 0: F1=0.322, sex AUROC=0.860, age MAE=8.888
- fold 1: F1=0.403, sex AUROC=0.951, age MAE=7.405
- fold 2: F1=0.442, sex AUROC=0.913, age MAE=7.546

### P4-fusion-balanced

_P4 saved scores + balanced logistic late fusion._

- fold 0: F1=0.302, sex AUROC=0.821, age MAE=9.916
- fold 1: F1=0.412, sex AUROC=0.960, age MAE=7.774
- fold 2: F1=0.410, sex AUROC=0.914, age MAE=8.169

### P4-fusion-logistic

_P4 saved scores + logistic late fusion._

- fold 0: F1=0.300, sex AUROC=0.821, age MAE=9.916
- fold 1: F1=0.399, sex AUROC=0.960, age MAE=7.774
- fold 2: F1=0.380, sex AUROC=0.914, age MAE=8.169

### P4-fusion-sgd

_P4 saved scores + SGD late fusion._

- fold 0: F1=0.240, sex AUROC=0.821, age MAE=9.916
- fold 1: F1=0.382, sex AUROC=0.960, age MAE=7.774
- fold 2: F1=0.302, sex AUROC=0.914, age MAE=8.169

### P4-pooling-mean

_Like P2 loss weights but **mean/mean** CpG→region and region→gene pooling; late-fusion test._

- fold 0: F1=0.300, sex AUROC=0.821, age MAE=9.916
- fold 1: F1=0.399, sex AUROC=0.960, age MAE=7.774
- fold 2: F1=0.380, sex AUROC=0.914, age MAE=8.169

### P5-epochs-30

_P2 loss weights, **max/max** pooling, **30-epoch** ceiling + early stop on validation tissue F1; late-fusion test._

- fold 0: F1=0.305, sex AUROC=0.820, age MAE=9.778
- fold 1: F1=0.404, sex AUROC=0.923, age MAE=7.838
- fold 2: F1=0.361, sex AUROC=0.906, age MAE=8.606

### P5-mean-epochs-30

_P5 schedule with **mean/mean** pooling._

- fold 0: F1=0.301, sex AUROC=0.815, age MAE=9.963
- fold 1: F1=0.400, sex AUROC=0.958, age MAE=7.666
- fold 2: F1=0.383, sex AUROC=0.907, age MAE=8.163

## Diagnosis

**Primary:** `task_competition`

P0 baseline tissue macro-F1 = 0.093 (7G replay).
P2 reweighted training (0.376) lifts tissue F1 >=0.05 vs P0 -> **task competition**.

## Checkpoint audit (trained cascade arms)

| Arm | Fold | Best epoch | Epochs run | Best val tissue F1 | Selection | Early stopped |
|-----|------|------------|------------|---------------------|-----------|---------------|
| P0-baseline | 0 | None | None | — | None | False |
| P0-baseline | 1 | None | None | — | None | False |
| P0-baseline | 2 | None | None | — | None | False |
| P2-end2end-tissue-weight | 0 | 15 | None | — | validation_tissue_macro_f1_then_age_mae | False |
| P2-end2end-tissue-weight | 1 | 12 | None | — | validation_tissue_macro_f1_then_age_mae | False |
| P2-end2end-tissue-weight | 2 | 9 | None | — | validation_tissue_macro_f1_then_age_mae | False |
| P4-pooling-mean | 0 | 9 | None | — | validation_tissue_macro_f1_then_age_mae | False |
| P4-pooling-mean | 1 | 12 | None | — | validation_tissue_macro_f1_then_age_mae | False |
| P4-pooling-mean | 2 | 12 | None | — | validation_tissue_macro_f1_then_age_mae | False |
| P5-epochs-30 | 0 | 9 | 17 | 0.231 | validation_tissue_macro_f1_then_age_mae | True |
| P5-epochs-30 | 1 | 4 | 12 | 0.116 | validation_tissue_macro_f1_then_age_mae | True |
| P5-epochs-30 | 2 | 4 | 12 | 0.241 | validation_tissue_macro_f1_then_age_mae | True |

## Study composition (tissue-labeled)

- fold 0 train: 5023 samples, 151 studies; top tissues: {46: 1540, 34: 497, 7: 345, 23: 345, 6: 275, 30: 152, 35: 140, 15: 109}
- fold 0 test: 2270 samples, 90 studies; top tissues: {46: 430, 22: 301, 24: 234, 7: 233, 35: 160, 15: 141, 34: 129, 26: 84}
- fold 1 train: 4164 samples, 159 studies; top tissues: {46: 926, 34: 421, 22: 330, 35: 261, 24: 246, 15: 205, 30: 170, 7: 167}
- fold 1 test: 3273 samples, 86 studies; top tissues: {46: 1072, 23: 365, 7: 328, 34: 208, 6: 201, 15: 95, 37: 90, 18: 90}
- fold 2 train: 4384 samples, 149 studies; top tissues: {46: 904, 7: 559, 22: 333, 23: 295, 15: 236, 6: 217, 34: 208, 35: 199}
- fold 2 test: 2323 samples, 86 studies; top tissues: {46: 568, 34: 292, 30: 230, 43: 109, 35: 101, 28: 99, 16: 84, 32: 82}

## Milestone 7 recommendation

**Proceed with narrowed hyperparameter grid** (P4-P5 + fusion solver sweep) before Milestone 7 OOF; cascade may be salvageable for tissue with fusion/loss fixes.

## Methodological note

P2 trains on MBS only but current test metrics late-fuse `[orphan_rbs | mbs | direct_contrib]`. P2 F1 above `C-mvalue-enet` does **not** prove MBS-only aggregation wins. Corrected benchmark: **7G′ Stage A** (`P2-G` … `C-mvalue-enet-G` on identical gene-linked CpGs).

## Product scores versus phenotype comparator

The same phenotype-trained model exports MBS (+ optional orphan/direct) **and** supports phenotype heads. Product export: 7F cascade topology. `direct_contrib.zarr` is diagnostic only; association needs `direct_cpg.zarr` (7G′ Stage B).

## Artifacts

- Config: `stage0_7g_cascade_tissue_probe`
- `arm_means.json`, `per_arm/*.json`, `figures/tissue_f1_bars.png`

Phase-2 is complete only when P4 and P5 each have all three fold artifacts. Do not infer their performance from P0-P3.
