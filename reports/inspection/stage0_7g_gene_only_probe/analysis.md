# 7G′ Stage A — gene-only MBS architecture selection

Primary metric: **`mbs_e2e`** tissue macro-F1 (end-to-end MBS heads; not late fusion).
Classical comparator: **`C-mvalue-*-G`** on identical `gene_cols`.

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
| `P2-G` | CascadeDeepSet on **gene-linked CpGs only** (`gene_cols`); max/max pooling; 15 epochs; P2 loss weights; primary metric **`mbs_e2e`** (end-to-end MBS heads). |
| `P4-G` | Like P2-G but **mean/mean** pooling; same gene-only panel and **`mbs_e2e`** metric. |
| `P5-G-max` | Gene-only cascade; max/max pooling; 30-epoch cap + early stop; **`mbs_e2e`** metric. |
| `P5-G-mean` | Gene-only cascade; mean/mean pooling; 30-epoch cap + early stop; **`mbs_e2e`** metric. |
| `C-mvalue-ridge-G` | Same as `C-mvalue-ridge` but only **`gene_cols`** (56k matched panel). |
| `C-mvalue-enet-G` | Elastic-net on **identical gene-linked CpGs** as neural `-G` arms (fair Stage A comparator). |
| `C-mvalue-hgb-G` | HGB on **`gene_cols`** only. |
| `C-mvalue-sva-G` | PCA-SVA + ridge on **`gene_cols`** only. |
| `C-mvalue-classical-G` | Bundle runner for all **`C-mvalue-*-G`** arms below on the same `gene_cols` panel. |
| `P2-orphan-ablation` | Full assignment (not gene-only): same P2 loss weights; compares **`fusion_full`** (orphan RBS + MBS + direct) vs **`fusion_mbs_direct`** (MBS + direct, orphan block dropped). |
| `C-mvalue-enetS` | **Stage B:** fold-safe stability-selected sparse CpG panel (outer-train only) → refit enet; one test eval per fold. |
| `N-cascade-S` | **Stage B:** locked Stage-A cascade hyperparameters on the **fold-selected** CpG panel (same loci as `C-mvalue-enetS` per fold). |
| `N-light-type` | **FlatDeepSetRegion** (planned): per-CpG `[M-value, regulatory-type one-hot, observed]` → pool by gene → MBS; lighter than full RBS→gene cascade. |
| `N-mbs-posthoc-full-fusion` | **Stage B post-hoc fusion:** MBS encoder trained once; CPU late fusion on **`[orphan_rbs | mbs | direct_contrib]`** (`fusion_full`). Not joint end-to-end. |
| `N-mbs-posthoc-mbs-direct` | **Stage B orphan ablation:** same encoder; CPU fusion on **`[mbs | direct_contrib]`** only (`fusion_mbs_direct`). |

### Evaluation modes (cascade metrics)

- **`mbs_e2e`**: End-to-end **`MultitaskHeads`** on MBS only (no orphan/direct columns) — **Stage A primary metric**.
- **`mbs_linear_probe`**: CPU linear probe fit on **saved MBS matrix** only (representation check, not fusion).
- **`mbs_enet`**: CPU elastic-net heads on the **same saved MBS** (age + tissue + sex). Same CascadeDeepSet weights as `mbs_e2e`; no encoder retrain.
- **`fusion_full`**: Late fusion on **`[orphan_rbs | mbs | direct_contrib]`** columns.
- **`fusion_mbs_direct`**: Late fusion on **`[mbs | direct_contrib]`** — orphan RBS ablation.

## Comparable ranking (panel × eval mode)

Compare **only within the same row group** (same panel and eval mode). Stage A primary metric is **`mbs_e2e`** on the **gene-linked** panel (test split only). 7G tissue probe P0–P5 used **late fusion (`fusion_full`)** on the **65k prefix**. Rows marked *invalid* used pre-fix `mbs_e2e` that scored train+validation+test together.

| Arm | Panel | Eval mode | Tissue macro-F1 | folds | Notes |
|-----|-------|-----------|----------------:|------:|-------|
| `P0-baseline` | 65k prefix | `fusion_full` | 0.093 (±0.014) | 3 |  |
| `P0-baseline` | 65k prefix | `mbs_e2e` | — | 0 | mbs_e2e missing eval_split=test (train+test leak) |
| `P2-end2end-tissue-weight` | 65k prefix | `fusion_full` | 0.376 (±0.059) | 3 |  |
| `P2-end2end-tissue-weight` | 65k prefix | `mbs_e2e` | — | 0 | mbs_e2e missing eval_split=test (train+test leak) |
| `P4-pooling-mean` | 65k prefix | `fusion_full` | 0.360 (±0.053) | 3 |  |
| `P4-pooling-mean` | 65k prefix | `mbs_e2e` | — | 0 | mbs_e2e missing eval_split=test (train+test leak) |
| `P5-epochs-30` | 65k prefix | `fusion_full` | 0.356 (±0.050) | 3 |  |
| `P5-epochs-30` | 65k prefix | `mbs_e2e` | — | 0 | mbs_e2e missing eval_split=test (train+test leak) |
| `P2-G` | gene-linked | `fusion_full` | 0.373 (±0.052) | 3 |  |
| `P2-G` | gene-linked | `mbs_e2e` | 0.373 (±0.038) | 3 |  |
| `P2-G` | gene-linked | `mbs_enet` | 0.385 (±0.053) | 3 |  |
| `P4-G` | gene-linked | `fusion_full` | 0.379 (±0.055) | 3 |  |
| `P4-G` | gene-linked | `mbs_e2e` | 0.370 (±0.059) | 3 |  |
| `P4-G` | gene-linked | `mbs_enet` | 0.379 (±0.059) | 3 |  |
| `P5-G-max` | gene-linked | `mbs_e2e` | 0.356 (±0.042) | 3 |  |
| `C-mvalue-enet` | 65k prefix | `classical` | 0.395 (±0.018) | 3 |  |
| `C-mvalue-enet-G` | gene-linked | `classical` | 0.388 (±0.018) | 3 |  |

**Fair pairs (examples):**

- **Late fusion, 65k:** `P2-end2end` vs `P4-pooling-mean` vs `C-mvalue-enet`.
- **Late fusion, gene-linked:** `P2-G` / `P4-G` `fusion_full` on `explicit_only`.
- **MBS e2e, gene-linked:** `P2-G` vs `P4-G` vs `P5-G-max` vs `C-mvalue-enet-G` (Stage A lock).
- **MBS enet readout:** same frozen MBS as e2e, elastic-net heads (not a lock substitute).

## Task comparison (tissue / age / sex)

Same **`explicit_only`** gene-linked panel and outer **test** folds. Compare rows as alternative **readouts** of one encoder (`mbs_e2e` / `mbs_linear_probe` / `mbs_enet`) versus classical models on raw CpG M-values. `mbs_e2e` sex AUROC is unavailable (class argmax only). Classical enet age is blanked when SGD MAE exploded. Horvath-style clocks are not in this table.

| Arm | Readout | Tissue F1 | Age MAE | Age R² | Sex AUROC | Sex F1 | folds |
|-----|---------|----------:|--------:|-------:|----------:|-------:|------:|
| `C-mvalue-enet-G` | `classical` | 0.388 (±0.018) | — | — | 0.882 (±0.033) | — | 3 |
| `P2-G` | `mbs_enet` | 0.385 (±0.053) | 14.393 (±0.655) | 0.452 (±0.054) | 0.765 (±0.027) | 0.702 (±0.026) | 3 |
| `P4-G` | `mbs_linear_probe` | 0.379 (±0.055) | 11.609 (±1.405) | 0.611 (±0.067) | 0.800 (±0.055) | 0.732 (±0.046) | 3 |
| `P4-G` | `mbs_enet` | 0.379 (±0.059) | 16.011 (±3.018) | 0.350 (±0.141) | 0.789 (±0.057) | 0.723 (±0.049) | 3 |
| `P2-G` | `mbs_e2e` | 0.373 (±0.038) | 15.637 (±0.787) | 0.315 (±0.140) | — | 0.708 (±0.030) | 3 |
| `P2-G` | `mbs_linear_probe` | 0.373 (±0.052) | 13.393 (±1.104) | 0.488 (±0.110) | 0.784 (±0.020) | 0.717 (±0.021) | 3 |
| `P5-G-max` | `mbs_linear_probe` | 0.371 (±0.039) | 13.274 (±0.893) | 0.476 (±0.091) | 0.714 (±0.009) | 0.659 (±0.003) | 3 |
| `P4-G` | `mbs_e2e` | 0.370 (±0.059) | 20.380 (±2.156) | -0.036 (±0.059) | — | 0.664 (±0.038) | 3 |
| `P5-G-max` | `mbs_e2e` | 0.356 (±0.042) | 21.402 (±1.786) | -0.150 (±0.078) | — | 0.598 (±0.006) | 3 |
| `C-mvalue-sva-G` | `classical` | 0.348 (±0.028) | 12.920 (±4.941) | 0.083 (±0.814) | 0.851 (±0.071) | — | 3 |
| `C-mvalue-ridge-G` | `classical` | 0.337 (±0.040) | 6.489 (±0.907) | 0.856 (±0.033) | 0.904 (±0.056) | — | 3 |
| `C-mvalue-hgb-G` | `classical` | 0.114 (±0.081) | 9.066 (±1.022) | 0.753 (±0.052) | 0.938 (±0.059) | — | 3 |

**Readouts:** `mbs_e2e` = jointly trained neural heads on MBS (Stage A lock metric); `mbs_linear_probe` / `mbs_enet` = new sklearn heads on the **same frozen MBS**; `classical` = sklearn on gene-linked CpG M-values (no encoder).

## Cascade arms (gene-linked CpGs only)

Primary **`mbs_e2e`** (test split only); **`mbs_linear_probe`** and **`mbs_enet`** are readouts of the **same frozen MBS**. Contaminated pre-fix **`mbs_e2e`** shown as *invalid*.

| Arm | mbs_e2e F1 | linear probe F1 | mbs_enet F1 | age MAE (e2e) | sex AUROC (probe) | folds |
|-----|-----------:|----------------:|------------:|--------------:|------------------:|------:|
| P2-G | 0.373 (±0.038) | 0.373 | 0.385 (±0.053) | 15.637 | 0.784 | 3 |
| P4-G | 0.370 (±0.059) | 0.379 | 0.379 (±0.059) | 20.380 | 0.800 | 3 |
| P5-G-max | 0.356 (±0.042) | 0.371 | — | 21.402 | 0.714 | 3 |

## Classical arms (-G panel)

Same **51,375 gene-linked CpGs** as neural arms (`explicit_only`): ridge, elastic-net, HGB, PCA-SVA+ridge.

| Arm | tissue F1 | age MAE | age R² | sex AUROC |
|-----|----------:|--------:|-------:|----------:|
| C-mvalue-enet-G | 0.388 (±0.018) | — | — | 0.882 |
| C-mvalue-sva-G | 0.348 (±0.028) | 12.920 | 0.083 | 0.851 |
| C-mvalue-ridge-G | 0.337 (±0.040) | 6.489 | 0.856 | 0.904 |
| C-mvalue-hgb-G | 0.114 (±0.081) | 9.066 | 0.753 | 0.938 |

## Locked architecture (Stage B input)

- **Cascade arm:** `P2-G`
- **Pooling (CpG / region):** `max` / `max`
- **Epoch ceiling:** 15
- **Best classical (-G):** `C-mvalue-enet-G`
- **Cascade clearly ahead (≥0.03 tissue F1):** False

**Encoder parity recommended:** re-run **FlatDeepSet** and **HierarchicalDeepSet** on the same `gene_cols` before committing to cascade for Stage B / Milestone 7.

## Orphan RBS ablation (P2-orphan-ablation)

Compare **`fusion_full`** (orphan RBS + MBS + direct) vs **`fusion_mbs_direct`** (MBS + direct only).

| Mode | Mean tissue macro-F1 |
|------|---------------------:|
| fusion_full | 0.368 |
| fusion_mbs_direct | 0.368 |
| Δ (full − mbs_direct) | 0.000 |

**Orphan RBS effect is negligible** at this budget (|Δ| ≤ 0.01); Stage B should still report both fusion modes.

## Parallel / follow-on work

- **Stage A required GPU arms** (`P2-G`, `P4-G`, `P5-G-max`, `C-mvalue-*-G`) are complete on `explicit_only`. Optional `P5-G-mean` was not run.
- **Encoder parity (optional):** FlatDeepSet + HierarchicalDeepSet on same `gene_cols` if cascade does not lead classical by ≥0.03 F1.
- **Stage B (after lock):** fold-safe `C-mvalue-enetS`, `N-cascade-S`, `N-light-type`, `direct_cpg.zarr`, full-model fusion arms.

## Next

- Stage B: fold-safe `C-mvalue-enetS`, `N-cascade-S`, `N-light-type` (FlatDeepSetRegion), `N-mbs-posthoc-full-fusion` / `N-mbs-posthoc-mbs-direct`, plus `direct_cpg.zarr`.
- Milestone **7** 5×6 OOF remains blocked until Stage B completes.
