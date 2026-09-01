# Stage 0 - 7G cascade tissue probe (P0-P3)

Frozen split `hub-ats-7e-3fold-v1`; budget **65536 loci / 15 epochs / 1 restart**.

**Phase-2 status:** P4/P5 code and configs are implemented, but this committed
report contains **P0–P3 results only**. P4/P5 performance is unknown until each
has three completed fold artifacts on the data/GPU host.

## Per-arm summary

| Arm | Tissue macro-F1 | Balanced acc | Sex AUROC | Age MAE | Age R² |
|-----|-----------------|--------------|-----------|---------|--------|
| P0-baseline | 0.093 | 0.125 | 0.911 | 8.441 | 0.782 |
| P1-fusion-tissue-heavy | 0.124 | 0.131 | 0.911 | 8.441 | 0.782 |
| P2-end2end-tissue-weight | 0.376 | 0.422 | 0.885 | 8.682 | 0.767 |
| P3-region-head-bypass | 0.389 | 0.424 | 0.908 | 7.946 | 0.802 |

## Per-fold tissue macro-F1

### P0-baseline
- fold 0: F1=0.076, sex AUROC=0.806, age MAE=9.653
- fold 1: F1=0.104, sex AUROC=0.982, age MAE=7.566
- fold 2: F1=0.098, sex AUROC=0.943, age MAE=8.105

### P1-fusion-tissue-heavy
- fold 0: F1=0.119, sex AUROC=0.806, age MAE=9.653
- fold 1: F1=0.119, sex AUROC=0.982, age MAE=7.566
- fold 2: F1=0.134, sex AUROC=0.943, age MAE=8.105

### P2-end2end-tissue-weight
- fold 0: F1=0.308, sex AUROC=0.818, age MAE=9.923
- fold 1: F1=0.407, sex AUROC=0.932, age MAE=7.828
- fold 2: F1=0.412, sex AUROC=0.907, age MAE=8.295

### P3-region-head-bypass
- fold 0: F1=0.322, sex AUROC=0.860, age MAE=8.888
- fold 1: F1=0.403, sex AUROC=0.951, age MAE=7.405
- fold 2: F1=0.442, sex AUROC=0.913, age MAE=7.546

## Paired comparison with locked 7G tissue comparator

| Fold | C-mvalue-enet | P2 cascade | P2 − C | P3 region mean |
|------|---------------|------------|--------|----------------|
| 0 | 0.303 | 0.308 | +0.005 | 0.322 |
| 1 | 0.335 | 0.407 | +0.072 | 0.403 |
| 2 | 0.364 | 0.412 | +0.048 | 0.442 |
| Mean | 0.334 | 0.376 | +0.042 | 0.389 |

This is encouraging, but it is three folds, one restart and a post-7G targeted
loss adjustment. It supports “salvageable,” not a new locked winner.

## Diagnosis

**Primary:** `task_competition`

P0 baseline tissue macro-F1 = 0.093 (7G replay).
P2 reweighted training (0.376) lifts tissue F1 >=0.05 vs P0 -> **task competition**.

## Checkpoint audit (P0 / P2)

| Arm | Fold | Best epoch | Selection |
|-----|------|------------|-----------|
| P0-baseline | 0 | None | None |
| P0-baseline | 1 | None | None |
| P0-baseline | 2 | None | None |
| P2-end2end-tissue-weight | 0 | 15 | validation_tissue_macro_f1_then_age_mae |
| P2-end2end-tissue-weight | 1 | 12 | validation_tissue_macro_f1_then_age_mae |
| P2-end2end-tissue-weight | 2 | 9 | validation_tissue_macro_f1_then_age_mae |

## Study composition (tissue-labeled)

- fold 0 train: 5023 samples, 151 studies; top tissues: {46: 1540, 34: 497, 7: 345, 23: 345, 6: 275, 30: 152, 35: 140, 15: 109}
- fold 0 test: 2270 samples, 90 studies; top tissues: {46: 430, 22: 301, 24: 234, 7: 233, 35: 160, 15: 141, 34: 129, 26: 84}
- fold 1 train: 4164 samples, 159 studies; top tissues: {46: 926, 34: 421, 22: 330, 35: 261, 24: 246, 15: 205, 30: 170, 7: 167}
- fold 1 test: 3273 samples, 86 studies; top tissues: {46: 1072, 23: 365, 7: 328, 34: 208, 6: 201, 15: 95, 37: 90, 18: 90}
- fold 2 train: 4384 samples, 149 studies; top tissues: {46: 904, 7: 559, 22: 333, 23: 295, 15: 236, 6: 217, 34: 208, 35: 199}
- fold 2 test: 2323 samples, 86 studies; top tissues: {46: 568, 34: 292, 30: 230, 43: 109, 35: 101, 28: 99, 16: 84, 32: 82}

## Milestone 7 recommendation

**Proceed with narrowed hyperparameter grid** (P4-P5 + fusion solver sweep) before Milestone 7 OOF; cascade may be salvageable for tissue with fusion/loss fixes.

## Product scores versus phenotype comparator

Product association export remains the deepMAT cascade: sample×gene MBS,
qualified orphan RBS kept one region per column, and indexed direct CpGs. The
current `direct_contrib.zarr` is a phenotype diagnostic, not yet the required
sample×direct-CpG association block. Tissue benchmarking is separate:
`C-mvalue-enet` remains the locked classical comparator (7G F1 0.334) until
the same-panel 7H comparison. See ADR 0010.

Current P2 F1 0.376 is numerically above the original C-mvalue-enet mean, but it
was obtained after a targeted tissue-loss investigation. It does not replace
the locked comparator or establish a final winner without P4/P5 and the
fold-selected same-panel benchmark.

## Artifacts

- Config: `stage0_7g_cascade_tissue_probe`
- `arm_means.json`, `per_arm/*.json`, `figures/tissue_f1_bars.png`

Phase-2 command after the implementation commit:

```bash
bash scripts/run_7g_cascade_tissue_probe.sh --device cuda \
  --arm P2-fusion-balanced \
  --arm P4-pooling-mean --arm P4-fusion-balanced \
  --arm P5-epochs-30 --arm P5-fusion-balanced
```

Do not mark Phase 2 complete or lock OOF hyperparameters until the regenerated
report contains all P4/P5 folds.
