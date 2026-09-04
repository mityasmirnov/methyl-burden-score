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
| `N-light-gene-mean` | One-hop annotated FlatDeepSetRegion with gene **mean** pooling. |
| `P2-G` | CascadeDeepSet on **gene-linked CpGs only** (`gene_cols`); max/max pooling; 15 epochs; P2 loss weights; primary metric **`mbs_e2e`** (end-to-end MBS heads). |
| `P4-G` | Like P2-G but **mean/mean** pooling; same gene-only panel and **`mbs_e2e`** metric. |
| `P5-G-max` | Gene-only cascade; max/max pooling; 30-epoch cap + early stop; **`mbs_e2e`** metric. **Inactive** for Stage A DeepRVAT screen (historical only). |
| `N-cascade-scalar-mean-max` | Scalar RBS cascade; **mean** CpG→region + **max** region→gene; Stage A screen. |
| `N-cascade-scalar-max-mean` | Scalar RBS cascade; **max** CpG→region + **mean** region→gene; Stage A screen. |
| `N-cascade-vector-mean-max` | Vector-region cascade: pool region **embeddings** (mean CpG→region, max region→gene) then ``rho_G`` → MBS. |
| `N-cascade-vector-max-max` | Vector-region cascade with max/max pooling of latent region embeddings. |
| `N-light-gene-max` | One-hop annotated FlatDeepSetRegion: `[M, gene-role, CGI context, regulatory, flags]` → gene **max** → MBS. |
| `C-mvalue-ridge-G` | Same as `C-mvalue-ridge` but only **`gene_cols`** (56k matched panel). |
| `C-mvalue-enet-G` | Elastic-net on **identical gene-linked CpGs** as neural `-G` arms (fair Stage A comparator). |
| `C-mvalue-hgb-G` | HGB on **`gene_cols`** only. |
| `C-mvalue-sva-G` | PCA-SVA + ridge on **`gene_cols`** only. |
| `C-mvalue-classical-G` | Bundle runner for all **`C-mvalue-*-G`** arms below on the same `gene_cols` panel. |
| `P2-orphan-ablation` | Full assignment (not gene-only): same P2 loss weights; compares **`fusion_full`** (orphan RBS + MBS + direct) vs **`fusion_mbs_direct`** (MBS + direct, orphan block dropped). |
| `N-cascade-scalar-max-max` | Alias of **`P2-G`**: scalar RBS cascade, CpG max → gene max. |
| `N-cascade-scalar-mean-mean` | Alias of **`P4-G`**: scalar RBS cascade, CpG mean → gene mean. |
| `rbs_linear_probe` | See docs/ARM_GLOSSARY.md (no inline entry for `rbs_linear_probe`). |
| `rbs_enet` | See docs/ARM_GLOSSARY.md (no inline entry for `rbs_enet`). |
| `C-mvalue-enetS` | **Stage B:** fold-safe stability-selected sparse CpG panel (outer-train only) → refit enet; one test eval per fold. |
| `N-cascade-S` | **Stage B:** locked Stage-A cascade hyperparameters on the **fold-selected** CpG panel (same loci as `C-mvalue-enetS` per fold). |
| `N-light-type` | Legacy Stage B name for FlatDeepSetRegion; prefer **`N-light-gene-max`** / **`N-light-gene-mean`** in Stage A screen. |
| `N-mbs-posthoc-full-fusion` | **Stage B post-hoc fusion:** MBS encoder trained once; CPU late fusion on **`[orphan_rbs | mbs | direct_contrib]`** (`fusion_full`). Not joint end-to-end. |
| `N-mbs-posthoc-mbs-direct` | **Stage B orphan ablation:** same encoder; CPU fusion on **`[mbs | direct_contrib]`** only (`fusion_mbs_direct`). |

### Evaluation modes (cascade metrics)

- **`mbs_e2e`**: End-to-end **`MultitaskHeads`** on MBS only (no orphan/direct columns) — **Stage A primary metric**.
- **`mbs_linear_probe`**: CPU linear probe fit on **saved MBS matrix** only (representation check, not fusion).
- **`mbs_enet`**: CPU elastic-net heads on the **same saved MBS** (age + tissue + sex). Same CascadeDeepSet weights as `mbs_e2e`; no encoder retrain.
- **`rbs_linear_probe`**: CPU linear probe on frozen **gene-linked RBS** (`all_gene_rbs.zarr`); diagnoses loss before vs after gene pooling.
- **`rbs_enet`**: CPU elastic-net on frozen gene-linked RBS (same matrix as `rbs_linear_probe`).
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

Same **`explicit_only`** gene-linked panel and outer **test** folds. Compare rows as alternative **readouts** of one encoder (`mbs_e2e` / `mbs_linear_probe` / `mbs_enet` / `rbs_*`) versus classical models on raw CpG M-values. Prefer `mbs_e2e` sex AUROC when present (proba path); `rbs_*` diagnose loss before vs after gene pooling. Classical enet age uses Huber SGD elastic-net (year-scale target + eta0=1e-4); unscaled squared-error SGD exploded on this panel. Horvath-style clocks are not in this table.

| Arm | Readout | Tissue F1 | Age MAE | Age R² | Sex AUROC | Sex F1 | folds |
|-----|---------|----------:|--------:|-------:|----------:|-------:|------:|
| `N-light-gene-max` | `mbs_enet_nested` | 0.402 (±0.039) | 10.479 (±1.981) | 0.664 (±0.086) | 0.733 (±0.016) | 0.699 (±0.013) | 3 |
| `N-cascade-scalar-mean-max` | `mbs_enet` | 0.397 (±0.024) | 15.055 (±0.082) | 0.362 (±0.056) | 0.628 (±0.026) | 0.589 (±0.046) | 2 |
| `N-light-gene-mean` | `mbs_linear_probe` | 0.393 (±0.056) | 11.026 (±1.064) | 0.639 (±0.088) | 0.766 (±0.043) | 0.704 (±0.038) | 3 |
| `C-mvalue-enet-G` | `classical` | 0.388 (±0.018) | 8.150 (±0.714) | 0.797 (±0.039) | 0.882 (±0.033) | — | 3 |
| `N-light-gene-mean` | `mbs_enet_nested` | 0.387 (±0.043) | 10.320 (±1.282) | 0.675 (±0.089) | 0.746 (±0.005) | 0.714 (±0.023) | 3 |
| `P2-G` | `mbs_enet` | 0.385 (±0.053) | 14.393 (±0.655) | 0.452 (±0.054) | 0.765 (±0.027) | 0.702 (±0.026) | 3 |
| `N-light-gene-max` | `mbs_enet` | 0.382 (±0.056) | 16.469 (±2.219) | 0.293 (±0.143) | 0.736 (±0.055) | 0.673 (±0.064) | 3 |
| `N-light-gene-mean` | `mbs_enet` | 0.382 (±0.055) | 15.324 (±2.687) | 0.395 (±0.113) | 0.730 (±0.051) | 0.672 (±0.051) | 3 |
| `P4-G` | `mbs_linear_probe` | 0.379 (±0.055) | 11.609 (±1.405) | 0.611 (±0.067) | 0.800 (±0.055) | 0.732 (±0.046) | 3 |
| `P4-G` | `mbs_enet` | 0.379 (±0.059) | 16.011 (±3.018) | 0.350 (±0.141) | 0.789 (±0.057) | 0.723 (±0.049) | 3 |
| `N-light-gene-mean` | `mbs_e2e` | 0.378 (±0.060) | 17.095 (±0.985) | 0.212 (±0.134) | 0.655 (±0.057) | 0.601 (±0.047) | 3 |
| `N-light-gene-max` | `mbs_linear_probe` | 0.375 (±0.041) | 12.048 (±0.940) | 0.577 (±0.020) | 0.746 (±0.069) | 0.685 (±0.057) | 3 |
| `P2-G` | `mbs_e2e` | 0.373 (±0.038) | 15.637 (±0.787) | 0.315 (±0.140) | — | 0.708 (±0.030) | 3 |
| `N-cascade-scalar-mean-max` | `rbs_linear_probe` | 0.373 (±0.045) | 15.231 (±2.720) | 0.273 (±0.294) | 0.764 (±0.019) | 0.699 (±0.013) | 3 |
| `P2-G` | `mbs_linear_probe` | 0.373 (±0.052) | 13.393 (±1.104) | 0.488 (±0.110) | 0.784 (±0.020) | 0.717 (±0.021) | 3 |
| `P5-G-max` | `mbs_linear_probe` | 0.371 (±0.039) | 13.274 (±0.893) | 0.476 (±0.091) | 0.714 (±0.009) | 0.659 (±0.003) | 3 |
| `P4-G` | `mbs_e2e` | 0.370 (±0.059) | 20.380 (±2.156) | -0.036 (±0.059) | — | 0.664 (±0.038) | 3 |
| `N-cascade-vector-mean-max` | `rbs_linear_probe` | 0.368 (±0.044) | 10.459 (±0.402) | 0.686 (±0.026) | 0.842 (±0.054) | 0.774 (±0.041) | 3 |
| `N-cascade-vector-max-max` | `mbs_linear_probe` | 0.367 (±0.046) | 14.857 (±1.089) | 0.367 (±0.048) | 0.687 (±0.029) | 0.638 (±0.019) | 3 |
| `N-cascade-vector-max-max` | `mbs_enet` | 0.366 (±0.048) | 15.345 (±2.305) | 0.386 (±0.095) | 0.672 (±0.029) | 0.628 (±0.020) | 3 |
| `N-cascade-vector-mean-max` | `mbs_enet` | 0.362 (±0.062) | 16.388 (±2.597) | 0.328 (±0.102) | 0.664 (±0.047) | 0.621 (±0.049) | 3 |
| `N-cascade-vector-mean-max` | `mbs_linear_probe` | 0.360 (±0.040) | 15.051 (±0.542) | 0.349 (±0.061) | 0.697 (±0.049) | 0.644 (±0.031) | 3 |
| `P5-G-max` | `mbs_e2e` | 0.356 (±0.042) | 21.402 (±1.786) | -0.150 (±0.078) | — | 0.598 (±0.006) | 3 |
| `N-cascade-scalar-mean-max` | `mbs_linear_probe` | 0.352 (±0.033) | 14.996 (±1.818) | 0.324 (±0.211) | 0.695 (±0.028) | 0.643 (±0.020) | 3 |
| `N-cascade-scalar-mean-max` | `mbs_enet_nested` | 0.349 (±0.037) | 13.738 (±0.947) | 0.269 (±0.167) | 0.673 (±0.040) | 0.642 (±0.027) | 3 |
| `C-mvalue-sva-G` | `classical` | 0.348 (±0.028) | 12.920 (±4.941) | 0.083 (±0.814) | 0.851 (±0.071) | — | 3 |
| `N-cascade-scalar-mean-max` | `mbs_e2e` | 0.346 (±0.050) | 20.411 (±1.502) | -0.049 (±0.004) | 0.683 (±0.029) | 0.623 (±0.021) | 3 |
| `N-cascade-vector-max-max` | `mbs_e2e` | 0.343 (±0.063) | 21.458 (±2.454) | -0.153 (±0.105) | 0.671 (±0.033) | 0.584 (±0.027) | 3 |
| `N-cascade-vector-mean-max` | `mbs_e2e` | 0.337 (±0.036) | 22.753 (±4.259) | -0.233 (±0.215) | 0.665 (±0.048) | 0.601 (±0.032) | 3 |
| `C-mvalue-ridge-G` | `classical` | 0.337 (±0.040) | 6.489 (±0.907) | 0.856 (±0.033) | 0.904 (±0.056) | — | 3 |
| `N-light-gene-max` | `mbs_e2e` | 0.336 (±0.047) | 21.593 (±4.138) | -0.263 (±0.580) | 0.624 (±0.039) | 0.547 (±0.040) | 3 |
| `N-cascade-vector-max-max` | `rbs_linear_probe` | 0.329 (±0.074) | 12.102 (±2.226) | 0.592 (±0.095) | 0.803 (±0.111) | 0.734 (±0.097) | 3 |
| `N-cascade-vector-mean-max` | `rbs_enet` | 0.316 (±0.050) | 19.667 (±2.204) | 0.078 (±0.075) | 0.837 (±0.047) | 0.720 (±0.030) | 3 |
| `N-cascade-scalar-max-mean` | `mbs_linear_probe` | 0.309 (±0.000) | 13.338 (±0.000) | 0.554 (±0.000) | 0.750 (±0.000) | 0.677 (±0.000) | 1 |
| `N-cascade-scalar-max-mean` | `rbs_linear_probe` | 0.298 (±0.000) | 15.308 (±0.000) | 0.370 (±0.000) | 0.753 (±0.000) | 0.681 (±0.000) | 1 |
| `N-cascade-scalar-max-mean` | `mbs_e2e` | 0.294 (±0.000) | 22.231 (±0.000) | -0.117 (±0.000) | 0.679 (±0.000) | 0.609 (±0.000) | 1 |
| `N-cascade-vector-max-max` | `rbs_enet` | 0.228 (±0.123) | 20.158 (±2.365) | 0.024 (±0.086) | 0.800 (±0.106) | 0.666 (±0.129) | 3 |
| `C-mvalue-hgb-G` | `classical` | 0.114 (±0.081) | 9.066 (±1.022) | 0.753 (±0.052) | 0.938 (±0.059) | — | 3 |

**Readouts:** `mbs_e2e` = jointly trained neural heads on MBS (ATS screen primary); `mbs_linear_probe` / `mbs_enet` = new sklearn heads on the **same frozen MBS**; `rbs_linear_probe` / `rbs_enet` = frozen **gene-linked RBS** (pre–gene-pool); `classical` = sklearn on gene-linked CpG M-values (no encoder). **folds** = number of folds that actually contain that readout (±0.000 with folds=1 means a single fold, not three identical scores).

## Three-task Pareto (`mbs_e2e` + classical)

Non-dominated on tissue macro-F1 (↑), age MAE (↓), sex AUROC (↑). Do **not** pick a winner on tissue alone.

| Arm | Readout | Tissue F1 | Age MAE | Sex AUROC |
|-----|---------|----------:|--------:|----------:|
| `C-mvalue-enet-G` | `classical` | 0.388 | 8.150 | 0.882 |
| `C-mvalue-ridge-G` | `classical` | 0.337 | 6.489 | 0.904 |
| `C-mvalue-hgb-G` | `classical` | 0.114 | 9.066 | 0.938 |


## Architecture questions (Stage A screen)

1. **CpG → region pool (mean vs max):** `mean-max` tissue F1=0.346 vs `max-max` 0.373; age MAE 20.411 vs 15.637. Prefer **`P2-G`** on this matched slice (check Pareto).
2. **Region → gene pool (mean vs max):** `max-mean` tissue F1=0.294 vs `max-max` 0.373; age MAE 22.231 vs 15.637. Prefer **`P2-G`** on this matched slice (check Pareto).
3. **Does scalar RBS discard information?** Vector arm `0.337` tissue vs P2 `0.373`; if vector does not beat scalar on age/sex, bottleneck is elsewhere.
4. **Gene pooling vs RBS:** `N-cascade-vector-mean-max` `rbs_*` tissue F1=0.368, age MAE=10.459, sex AUROC=0.842; same-arm MBS probe tissue=0.362, age=16.388, sex=0.664. Gene pooling is near-neutral on tissue; **age/sex often better on RBS** (pre–gene-pool), so some phenotype signal is lost at region→gene. Classical enet age MAE=8.150 remains the age ceiling.
5. **One-hop vs cascade:** One-hop `N-light-gene-max` tissue=0.336 / age=21.593 vs P2-G 0.373 / 15.637.
6. **One-scalar-per-gene bottleneck:** Gene aggregation still trails classical on age/sex; one scalar MBS/gene is **not yet adequate** unless a screen arm closes the gap.
7. **Best performance/compute:** Prefer landed P2/P4 (15 ep) as the current ATS **reference**, not a pooling lock. Do **not** promote unmatched Tier-1 (5 ep) arms against 15-ep P2. Next gate is the **matched 16-epoch promotion screen**; age-primary seed-mask waits on those decision rules.

## Training epochs (ceiling / ran / best)

Ceiling is the configured `max_epochs` (Tier-1 screen note for N-light / mixed/vector arms). **Ran** is how many epochs the trainer completed (early stop may cut short). **Best** is the checkpoint selected by `validation_tissue_macro_f1_then_age_mae` (used for test `mbs_e2e`). Prefer actual best/ran over the ceiling label — do not stamp a hard-coded 5-epoch N-light label onto longer runs.

| Arm | Ceiling | Epochs ran (per fold) | Best epoch (per fold) | folds |
|-----|--------:|----------------------:|----------------------:|------:|
| `P2-G` | 15 | — | 15,9,6 (μ=10.0) | 3 |
| `P4-G` | 15 | — | 12,6,15 (μ=11.0) | 3 |
| `P5-G-max` | 30 | — | 13,9,13 (μ=11.7) | 3 |
| `N-cascade-scalar-mean-max` | 16 | 13,14,16 (μ=14.3) | 8,9,15 (μ=10.7) | 3 |
| `N-cascade-scalar-max-mean` | 16 | 16 | 13 | 1 |
| `N-cascade-vector-mean-max` | 16 | 5,5,5 (μ=5.0) | 5,5,4 (μ=4.7) | 3 |
| `N-cascade-vector-max-max` | 5 | 5,5,5 (μ=5.0) | 5,4,5 (μ=4.7) | 3 |
| `N-light-gene-max` | 16 | 21,14,16 (μ=17.0) | 16,9,16 (μ=13.7) | 3 |
| `N-light-gene-mean` | 16 | 16,16,16 (μ=16.0) | 16,14,15 (μ=15.0) | 3 |

## Cascade arms (gene-linked CpGs only)

Primary **`mbs_e2e`** (test split only); **`mbs_linear_probe`** and **`mbs_enet`** are readouts of the **same frozen MBS**; **`rbs_linear_probe`** / **`rbs_enet`** use gene-linked RBS (`all_gene_rbs.zarr`). Contaminated pre-fix **`mbs_e2e`** shown as *invalid*. **Best ep** = checkpoint epoch used for test eval; **ran** = epochs completed.

| Arm | mbs_e2e F1 | linear probe F1 | mbs_enet F1 | age MAE (e2e) | sex AUROC (probe) | best ep | ran | folds |
|-----|-----------:|----------------:|------------:|--------------:|------------------:|--------:|----:|------:|
| N-light-gene-mean | 0.378 (±0.060) | 0.393 | 0.382 (±0.055) | 17.095 | 0.766 | 16,14,15 (μ=15.0) | 16,16,16 (μ=16.0) | 3 |
| P2-G | 0.373 (±0.038) | 0.373 | 0.385 (±0.053) | 15.637 | 0.784 | 15,9,6 (μ=10.0) | — | 3 |
| P4-G | 0.370 (±0.059) | 0.379 | 0.379 (±0.059) | 20.380 | 0.800 | 12,6,15 (μ=11.0) | — | 3 |
| P5-G-max | 0.356 (±0.042) | 0.371 | — | 21.402 | 0.714 | 13,9,13 (μ=11.7) | — | 3 |
| N-cascade-scalar-mean-max | 0.346 (±0.050) | 0.352 | 0.397 (±0.024) [2/3] | 20.411 | 0.695 | 8,9,15 (μ=10.7) | 13,14,16 (μ=14.3) | 3 |
| N-cascade-vector-max-max | 0.343 (±0.063) | 0.367 | 0.366 (±0.048) | 21.458 | 0.687 | 5,4,5 (μ=4.7) | 5,5,5 (μ=5.0) | 3 |
| N-cascade-vector-mean-max | 0.337 (±0.036) | 0.360 | 0.362 (±0.062) | 22.753 | 0.697 | 5,5,4 (μ=4.7) | 5,5,5 (μ=5.0) | 3 |
| N-light-gene-max | 0.336 (±0.047) | 0.375 | 0.382 (±0.056) | 21.593 | 0.746 | 16,9,16 (μ=13.7) | 21,14,16 (μ=17.0) | 3 |
| N-cascade-scalar-max-mean | 0.294 (±0.000) | 0.309 | — | 22.231 | 0.750 | 13 | 16 | 1 |

### RBS frozen readouts (screen cascade — `rbs_enet` / `rbs_linear`)

| Arm | `rbs_enet` tissue F1 | `rbs_enet` age MAE | `rbs_enet` sex AUROC | `rbs_linear` tissue F1 | `rbs_linear` age MAE | folds |
|-----|---------------------:|-------------------:|---------------------:|-----------------------:|---------------------:|------:|
| `N-cascade-scalar-mean-max` | — | — | — | 0.373 (±0.045) | 15.231 (±2.720) | 3 |
| `N-cascade-scalar-max-mean` | — | — | — | 0.298 (±0.000) | 15.308 (±0.000) | 1 |
| `N-cascade-vector-mean-max` | 0.316 (±0.050) | 19.667 (±2.204) | 0.837 (±0.047) | 0.368 (±0.044) | 10.459 (±0.402) | 3 |
| `N-cascade-vector-max-max` | 0.228 (±0.123) | 20.158 (±2.365) | 0.800 (±0.106) | 0.329 (±0.074) | 12.102 (±2.226) | 3 |

`rbs_enet` via `scripts/eval_mbs_enet_from_scores.py --which rbs` on saved `all_gene_rbs.zarr` (13 212 regions; no encoder retrain). Fixed `alpha=0.1` / `l1_ratio=0.5` **without** train-fold standardization is **diagnostic only**. Scalar arms: enet ≈/≥ linear on tissue and improves age. Vector arms: age collapses under that fixed enet while sex stays nearly unchanged — that is an over-strong / unscaled sparse penalty, **not** evidence the vector RBS representation is weak. Prefer `rbs_linear_probe` (and nested `rbs_enet_nested` once available) for vector RBS. P2-G `rbs_enet` not run (folds 1–2 lack `all_gene_rbs.zarr`).


## Classical arms (-G panel)

Same **51,375 gene-linked CpGs** as neural arms (`explicit_only`): ridge, elastic-net, HGB, PCA-SVA+ridge.

| Arm | tissue F1 | age MAE | age R² | sex AUROC |
|-----|----------:|--------:|-------:|----------:|
| C-mvalue-enet-G | 0.388 (±0.018) | 8.150 | 0.797 | 0.882 |
| C-mvalue-sva-G | 0.348 (±0.028) | 12.920 | 0.083 | 0.851 |
| C-mvalue-ridge-G | 0.337 (±0.040) | 6.489 | 0.856 | 0.904 |
| C-mvalue-hgb-G | 0.114 (±0.081) | 9.066 | 0.753 | 0.938 |

## Screen status (no architecture lock)

The ATS gene-only screen is **evidence**, not an architecture decision. **No cascade topology is locked.** Fold-selected-panel Stage B is **not** the next gate.

- **Best landed ATS cascade row:** `N-light-gene-mean` (pooling `n/a` / `mean`; configured ceiling 16 ep — see Training epochs for actual best/ran)
- **Best classical (-G):** `C-mvalue-enet-G`
- **Cascade clearly ahead (≥0.03 tissue F1):** False
- **Architecture locked:** `False`

Caveats: not ≥0.03 ahead of classical; tissue-primary loss; `P2-G` is the **current reference**, not a pooling lock. Unmatched Tier-1 (5 ep) cells must not be compared to 15-ep P2 as if budgets matched — run the 16-epoch promotion screen first.

## Orphan RBS ablation (P2-orphan-ablation)

Compare **`fusion_full`** (orphan RBS + MBS + direct) vs **`fusion_mbs_direct`** (MBS + direct only).

| Mode | Mean tissue macro-F1 |
|------|---------------------:|
| fusion_full | 0.368 |
| fusion_mbs_direct | 0.368 |
| Δ (full − mbs_direct) | 0.000 |

**Orphan RBS effect is negligible** at this budget (|Δ| ≤ 0.01); Stage B should still report both fusion modes.

## Annotation ablation grid (A0–A4, N0–N3)

Fold 0, `mean` pooling, ≤8 epochs, **two seeds** (primary + `-s2`) pooled. Bootstrap 95% CIs over seed runs. Primary metric `mbs_e2e` tissue macro-F1; linear probe is the representation check.

**Payloads found:** 18/18 seed runs under `per_arm/N-light-gene-ablation-*.json`.

**Note:** A4 ≈ N2 ≈ N3 while regulatory channels are zero (cCRE/DHS/ChromHMM not on disk). **`m_only` should lead** if annotations add noise under this budget.

| Arm | Features | Tissue e2e [95% CI] | Linear F1 [95% CI] | Age MAE (e2e) [95% CI] | Sex AUROC (e2e) [95% CI] |
|-----|----------|--------------------:|-------------------:|-----------------------:|-------------------------:|
| A0 | M only | 0.276 [0.275–0.276] | 0.350 [0.350–0.350] | 20.091 [20.051–20.131] | 0.639 [0.638–0.639] |
| A1 | M + gene role | 0.107 [0.107–0.107] | 0.316 [0.316–0.316] | 21.022 [20.992–21.053] | 0.595 [0.595–0.595] |
| A2 | M + CpG context | 0.170 [0.170–0.170] | 0.318 [0.317–0.319] | 23.251 [23.209–23.292] | 0.596 [0.596–0.597] |
| A3 | M + role + context | 0.168 [0.165–0.171] | 0.319 [0.318–0.319] | 21.563 [21.545–21.581] | 0.596 [0.592–0.600] |
| A4/A7 | All (regulatory zero) | 0.174 [0.173–0.175] | 0.319 [0.319–0.319] | 21.559 [21.422–21.696] | 0.590 [0.582–0.598] |

### Negative controls

| Arm | Features | Tissue e2e [95% CI] | Linear F1 [95% CI] | Age MAE (e2e) [95% CI] | Sex AUROC (e2e) [95% CI] |
|-----|----------|--------------------:|-------------------:|-----------------------:|-------------------------:|
| N0 | Observed flag only | 0.007 [0.006–0.008] | 0.011 [0.010–0.011] | 21.171 [21.116–21.227] | 0.540 [0.532–0.549] |
| N1 | Annotations only (no M) | 0.009 [0.009–0.009] | 0.019 [0.018–0.020] | 21.628 [21.623–21.633] | 0.531 [0.531–0.532] |
| N2 | Reg. permuted | 0.173 [0.172–0.174] | 0.319 [0.318–0.320] | 21.269 [21.233–21.306] | 0.586 [0.585–0.588] |
| N3 | All-zero regulatory | 0.169 [0.165–0.173] | 0.318 [0.318–0.319] | 21.455 [21.315–21.596] | 0.589 [0.586–0.592] |

**Takeaway:** best e2e tissue = **`A0` (M only)** at 0.276. `m_only` leads; gene-role/context do not help under this fold-0 budget. Negatives `obs_only` / `anno_only` should be near chance.


### Representation diagnostics (fold 0 mean across seeds)

Computed post-hoc from saved `scores/mbs.npy` (+ `mbs_present.npy`), checkpoint `head_state`, and Pearson r of per-sample mean MBS vs mean M-value over the gene-linked CpG panel (`sample_mean_m_gene_panel.npy`). Saturation = fraction of present scores ≤0.05 or ≥0.95; const-score = fraction of genes with SD < 1e-4 across samples.

| Arm | Gene-score SD | Saturation frac | Const-score frac | Corr w/ mean-M | Head ‖w‖₂ | Best ep |
|-----|:-------------:|:---------------:|:----------------:|:--------------:|:---------:|:-------:|
| A0 | 0.149 | 0.000 | 0.000 | -0.011 | 353.288 | 8 |
| A1 | 0.104 | 0.000 | 0.000 | 0.005 | 353.218 | 6 |
| A2 | 0.145 | 0.000 | 0.000 | 0.029 | 353.242 | 7 |
| A3 | 0.189 | 0.000 | 0.000 | 0.007 | 353.292 | 8 |
| A4/A7 | 0.188 | 0.000 | 0.000 | 0.007 | 353.293 | 8 |
| N0 | 0.000 | 0.000 | 1.000 | — | 353.163 | 5 |
| N1 | 0.042 | 0.000 | 0.033 | -0.093 | 353.184 | 7 |
| N2 | 0.188 | 0.000 | 0.000 | 0.009 | 353.292 | 8 |
| N3 | 0.188 | 0.000 | 0.000 | 0.008 | 353.292 | 8 |

**Repr read:** A0 gene-score SD≈0.149 (non-collapsed encoder); N0 const-score≈1 (obs-only scores collapsed — control OK); corr(mean MBS, panel mean-M)≈0 across arms (gene MBS ≠ bulk methylation intensity); no score saturation (not stuck at 0/1).


## Interpretation

### N-light is not collapsed

N-light mean improves from age MAE ~23.08 end-to-end to ~11.55 with a refitted linear head. The encoder contains information; the native optimisation / readout is poor.

### Annotation graph is populated; network input is weak

Audit: 51,375 unique CpGs; 57,430 locus–gene edges; 2,646 genes; 5,718 multi-gene CpGs; zero `other_gene` edges; all CpG-context categories populated. The short raw-concatenation ablation does **not** show annotations are uninformative — it shows **raw concatenation** hurts a short, tissue-primary, mean-pooling run. Implementation weaknesses (`gather_flat_region_features` / `FlatDeepSetRegion`):

- M-value + six gene-role one-hots + seven context one-hots + six regulatory slots + flags → one 24-d input;
- all six regulatory channels currently zero;
- `observed` effectively always one (unobserved edges dropped before encode);
- `gene_role_present` effectively constant on the Stage A graph;
- raw unnormalized M mixed with 0/1 annotations;
- global max/mean pool lets promoter and body cancel;
- fold-fitted robust-z fitted in `loop.py` but unused by the flat-region path.

Preferred one-hop (document only; do not train in this gate): gated embeddings `h_i = φ_M(z_i) + α_R E_R(role) + α_C E_C(context)` with `α` near 0, then **role-stratified** pools (promoter-core / proximal / 5′ / body / 3′) → gene embedding → scalar MBS.

### Vector vs scalar

Fair five-epoch mean→max: scalar tissue F1 ~0.331 vs vector ~0.337. Do **not** compare five-epoch vector to fifteen-epoch P2. Vector RBS **linear** probe before gene pool: age MAE ~10.46, sex AUROC ~0.842 — CpG→region works; failure is later (elementwise max/mean, no typed output channel, one scalar MBS). Fixed **`rbs_enet` on vector mean→max collapses age** (MAE ~19.67) while sex holds (~0.837 vs linear 0.842). That pattern is characteristic of an overly strong or poorly calibrated sparse penalty on a distributed age signal — **not** evidence that the vector RBS representation is weak. Prefer `rbs_linear_probe` / nested enet for vector diagnostics. Raising LR is not the first response.

### RBS → MBS and typed pooling

Cascade already adds a region-type embedding before RBS; the scalar path then pools only scalars (type discarded) and the vector path can still mix roles. CPU ablation **R0–R5** (`reports/inspection/stage0_7g_gene_only_probe/typed_rbs_pooling/`) finds typed max/mean (R1–R3) improve age MAE by ~3–4 y vs untyped R0 on vector-mean-max RBS, but the within-gene **role shuffle control does not collapse** — so the gain is not yet proof of biological role identity (extra channels / capacity may explain it). This diagnostic does **not** decide Stage B. A neural typed aggregator remains a pooling follow-up only if typed arms beat R0 on age **and** the shuffle control collapses.

### Next real gate

**Matched 16-epoch promotion screen** (one-hop max/mean, scalar mixed pools, vector mean→max) before declaring pooling winners or starting age-primary seed-mask. `P2-G` is the current reference, **not** a pooling lock. Fold-selected-panel Stage B stays blocked.

## Parallel / follow-on work

- **ATS Stage A Tier-1 screen + annotation ablations:** complete. Freeze **`P2-G` as current reference, not a pooling lock.**
- **Matched 16-epoch promotion screen (current GPU gate):** **`N-light-gene-max`** tissue e2e **0.336** (below P2 — do not rerun); **`N-light-gene-mean`** tissue e2e **0.378** (**within ~0.03 of P2** — `one_hop_mean_near_p2` fired); **`N-cascade-scalar-mean-max` 16-ep** done 3/3 (mean e2e ≈0.346). **Now:** scalar **max→mean** 16-ep (fold 0 done; fold 1+ training), then vector mean→max ×3 — [`milestone-7g-prime-16ep-promotion.md`](../../../docs/plans/milestone-7g-prime-16ep-promotion.md).
- **Post-hoc CPU enet:** light max/mean 16-ep already have nested `mbs_enet`. Launching / filling fixed+nested **mbs/rbs** enet on `scalar-mean-max-16ep` and available `scalar-max-mean-16ep` folds (`scratch/logs/16ep_posthoc_enet.log`; unit `mbs-16ep-enet`). Fixed `rbs_enet` remains diagnostic-only.
- **CPU typed-RBS ablation (R0–R5):** done; shuffle did not collapse → neural typed aggregator **not** promoted.
- **Age-primary seed-mask screen:** fold-0 panel audit **green** (`ok_for_seed_mask_gpu`); CUDA blocked only on 16-ep unlock — [`milestone-7g-prime-age-seed-mask.md`](../../../docs/plans/milestone-7g-prime-age-seed-mask.md).
- **Atlas association catalog:** done (SQL 013); non-blocking.
- **Stage B CpG-panel GPU:** blocked until seed-mask screen + typed-RBS diagnostics.

## Next

Ordered ops:

1. **Finish remaining 16-ep cascade queue** (scalar max→mean folds 1–2 → vector mean→max ×3); do not kill GPU 0 jobs.
2. **Let post-hoc enet finish** then re-sync `per_arm/` + `write_7g_gene_only_probe_report.py` / `apply_7g_16ep_decision.py`.
3. **Age-primary seed-mask GPU** (`scripts/run_7g_prime_seed_mask.py --device cuda --reuse-panels`) — only after `promotion_decision.json` unlocks; fold 0, two seeds, K=256; do **not** launch Stage B.
4. **Stage B** fold-selected CpG panel only after seed-mask screen.
5. **Milestone 7** 5×6 OOF after Stage B + `direct_cpg.zarr`.
