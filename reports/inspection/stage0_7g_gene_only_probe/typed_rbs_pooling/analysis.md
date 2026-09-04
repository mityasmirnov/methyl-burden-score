# Typed RBS pooling ablation (R0–R5)

Source run: `stage0-7g-gene-probe-vector-mean-max`. Presence-aware role features; Ridge/LogReg heads; train-fold normalization. **Amended promotion gate** (not age MAE alone).

| Arm | Tissue F1 | Age MAE | Sex AUROC | folds |
|-----|----------:|--------:|----------:|------:|
| `R0_max` | 0.326 (±0.041) | 14.660 (±1.538) | 0.719 (±0.047) | 3 |
| `R0_mean` | 0.289 (±0.042) | 14.599 (±1.459) | 0.731 (±0.063) | 3 |
| `R1` | 0.362 (±0.044) | 11.050 (±0.447) | 0.828 (±0.059) | 3 |
| `R2` | 0.366 (±0.043) | 11.020 (±0.419) | 0.829 (±0.060) | 3 |
| `R3` | 0.370 (±0.042) | 10.798 (±0.393) | 0.835 (±0.051) | 3 |
| `R4` | 0.347 (±0.047) | 12.606 (±0.864) | 0.808 (±0.074) | 3 |
| `R5` | 0.368 (±0.044) | 10.459 (±0.402) | 0.842 (±0.054) | 3 |
| `R1_shuffle` | 0.365 (±0.038) | 11.138 (±0.633) | 0.830 (±0.058) | 3 |

## Promotion gate (amended)

- Best typed arm: **`R3`**
- Age MAE Δ vs R0_max: **3.862** y (need ≈≥1)
- Age R² Δ: **0.222** (need ≈≥0.05)
- Tissue F1 Δ: **0.044** (loss ≤0.03)
- Sex AUROC Δ: **0.116** (loss ≤0.03)
- Shuffle control age MAE: **11.137672337936886** vs typed **10.798153059508705** (must clearly deteriorate)
- **Promote neural typed aggregator:** `False` (shuffle Δ=0.340 y; need ≥1)

Presence-aware encoding: `x=(RBS-0.5)*present` plus role-present flags / n_regions. Shuffle reports fraction of genes altered and columns whose role changed (see fold logs).
