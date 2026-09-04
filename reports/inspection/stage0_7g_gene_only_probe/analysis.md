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
| `P5-G-max` | Gene-only cascade; max/max pooling; 30-epoch cap + early stop; **`mbs_e2e`** metric. **Inactive** for Stage A DeepRVAT screen (historical only). |
| `N-cascade-scalar-mean-max` | Scalar RBS cascade; **mean** CpG→region + **max** region→gene; Stage A screen. |
| `N-cascade-scalar-max-mean` | Scalar RBS cascade; **max** CpG→region + **mean** region→gene; Stage A screen. |
| `N-cascade-vector-mean-max` | Vector-region cascade: pool region **embeddings** (mean CpG→region, max region→gene) then ``rho_G`` → MBS. |
| `N-cascade-vector-max-max` | Vector-region cascade with max/max pooling of latent region embeddings. |
| `N-light-gene-max` | One-hop annotated FlatDeepSetRegion: `[M, gene-role, CGI context, regulatory, flags]` → gene **max** → MBS. |
| `N-light-gene-mean` | One-hop annotated FlatDeepSetRegion with gene **mean** pooling. |
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
| `C-mvalue-enet-G` | `classical` | 0.388 (±0.018) | 8.150 (±0.714) | 0.797 (±0.039) | 0.882 (±0.033) | — | 3 |
| `P2-G` | `mbs_enet` | 0.385 (±0.053) | 14.393 (±0.655) | 0.452 (±0.054) | 0.765 (±0.027) | 0.702 (±0.026) | 3 |
| `N-cascade-scalar-mean-max` | `rbs_linear_probe` | 0.385 (±0.042) | 14.668 (±1.041) | 0.362 (±0.130) | 0.773 (±0.013) | 0.700 (±0.020) | 3 |
| `N-cascade-scalar-max-mean` | `rbs_linear_probe` | 0.383 (±0.066) | 13.507 (±2.237) | 0.469 (±0.135) | 0.811 (±0.055) | 0.738 (±0.045) | 3 |
| `P4-G` | `mbs_linear_probe` | 0.379 (±0.055) | 11.609 (±1.405) | 0.611 (±0.067) | 0.800 (±0.055) | 0.732 (±0.046) | 3 |
| `P4-G` | `mbs_enet` | 0.379 (±0.059) | 16.011 (±3.018) | 0.350 (±0.141) | 0.789 (±0.057) | 0.723 (±0.049) | 3 |
| `N-cascade-scalar-mean-max` | `mbs_linear_probe` | 0.375 (±0.029) | 14.382 (±1.173) | 0.390 (±0.161) | 0.711 (±0.031) | 0.654 (±0.024) | 3 |
| `P2-G` | `mbs_e2e` | 0.373 (±0.038) | 15.637 (±0.787) | 0.315 (±0.140) | — | 0.708 (±0.030) | 3 |
| `P2-G` | `mbs_linear_probe` | 0.373 (±0.052) | 13.393 (±1.104) | 0.488 (±0.110) | 0.784 (±0.020) | 0.717 (±0.021) | 3 |
| `P5-G-max` | `mbs_linear_probe` | 0.371 (±0.039) | 13.274 (±0.893) | 0.476 (±0.091) | 0.714 (±0.009) | 0.659 (±0.003) | 3 |
| `P4-G` | `mbs_e2e` | 0.370 (±0.059) | 20.380 (±2.156) | -0.036 (±0.059) | — | 0.664 (±0.038) | 3 |
| `N-cascade-vector-mean-max` | `rbs_linear_probe` | 0.368 (±0.044) | 10.459 (±0.402) | 0.686 (±0.026) | 0.842 (±0.054) | 0.774 (±0.041) | 3 |
| `N-cascade-scalar-max-mean` | `mbs_linear_probe` | 0.367 (±0.051) | 12.247 (±1.077) | 0.567 (±0.012) | 0.800 (±0.046) | 0.726 (±0.044) | 3 |
| `N-cascade-vector-max-max` | `mbs_linear_probe` | 0.367 (±0.046) | 14.857 (±1.089) | 0.367 (±0.048) | 0.687 (±0.029) | 0.638 (±0.019) | 3 |
| `N-cascade-vector-mean-max` | `mbs_linear_probe` | 0.360 (±0.040) | 15.051 (±0.542) | 0.349 (±0.061) | 0.697 (±0.049) | 0.644 (±0.031) | 3 |
| `N-cascade-scalar-max-mean` | `mbs_e2e` | 0.359 (±0.075) | 21.356 (±2.603) | -0.135 (±0.147) | 0.701 (±0.041) | 0.641 (±0.041) | 3 |
| `P5-G-max` | `mbs_e2e` | 0.356 (±0.042) | 21.402 (±1.786) | -0.150 (±0.078) | — | 0.598 (±0.006) | 3 |
| `C-mvalue-sva-G` | `classical` | 0.348 (±0.028) | 12.920 (±4.941) | 0.083 (±0.814) | 0.851 (±0.071) | — | 3 |
| `N-light-gene-mean` | `mbs_linear_probe` | 0.345 (±0.053) | 11.547 (±1.854) | 0.608 (±0.079) | 0.770 (±0.084) | 0.706 (±0.073) | 3 |
| `N-cascade-vector-max-max` | `mbs_e2e` | 0.343 (±0.063) | 21.458 (±2.454) | -0.153 (±0.105) | 0.671 (±0.033) | 0.584 (±0.027) | 3 |
| `N-cascade-vector-mean-max` | `mbs_e2e` | 0.337 (±0.036) | 22.753 (±4.259) | -0.233 (±0.215) | 0.665 (±0.048) | 0.601 (±0.032) | 3 |
| `C-mvalue-ridge-G` | `classical` | 0.337 (±0.040) | 6.489 (±0.907) | 0.856 (±0.033) | 0.904 (±0.056) | — | 3 |
| `N-cascade-scalar-mean-max` | `mbs_e2e` | 0.331 (±0.020) | 20.713 (±1.893) | -0.083 (±0.141) | 0.667 (±0.032) | 0.597 (±0.024) | 3 |
| `N-cascade-vector-max-max` | `rbs_linear_probe` | 0.329 (±0.074) | 12.102 (±2.226) | 0.592 (±0.095) | 0.803 (±0.111) | 0.734 (±0.097) | 3 |
| `N-light-gene-max` | `mbs_linear_probe` | 0.295 (±0.122) | 13.438 (±2.818) | 0.483 (±0.171) | 0.731 (±0.107) | 0.678 (±0.089) | 3 |
| `N-light-gene-mean` | `mbs_enet` | 0.278 (±0.000) | 20.882 (±0.000) | -0.008 (±0.000) | 0.760 (±0.000) | 0.553 (±0.000) | 3 |
| `N-light-gene-max` | `mbs_enet` | 0.131 (±0.000) | 20.872 (±0.000) | -0.022 (±0.000) | 0.627 (±0.000) | 0.536 (±0.000) | 3 |
| `N-light-gene-mean` | `mbs_e2e` | 0.126 (±0.064) | 23.082 (±2.533) | -0.406 (±0.533) | 0.591 (±0.014) | 0.462 (±0.091) | 3 |
| `N-light-gene-max` | `mbs_e2e` | 0.122 (±0.102) | 21.888 (±0.553) | -0.232 (±0.120) | 0.575 (±0.022) | 0.458 (±0.117) | 3 |
| `C-mvalue-hgb-G` | `classical` | 0.114 (±0.081) | 9.066 (±1.022) | 0.753 (±0.052) | 0.938 (±0.059) | — | 3 |

**Readouts:** `mbs_e2e` = jointly trained neural heads on MBS (Stage A lock metric); `mbs_linear_probe` / `mbs_enet` = new sklearn heads on the **same frozen MBS**; `rbs_linear_probe` / `rbs_enet` = frozen **gene-linked RBS** (pre–gene-pool); `classical` = sklearn on gene-linked CpG M-values (no encoder).

## Three-task Pareto (`mbs_e2e` + classical)

Non-dominated on tissue macro-F1 (↑), age MAE (↓), sex AUROC (↑). Do **not** pick a winner on tissue alone.

| Arm | Readout | Tissue F1 | Age MAE | Sex AUROC |
|-----|---------|----------:|--------:|----------:|
| `C-mvalue-enet-G` | `classical` | 0.388 | 8.150 | 0.882 |
| `C-mvalue-ridge-G` | `classical` | 0.337 | 6.489 | 0.904 |
| `C-mvalue-hgb-G` | `classical` | 0.114 | 9.066 | 0.938 |


## Interpretation (2026-09-04)

Tier-1 DeepRVAT screen + fold-0 annotation ablations are in. **No new arm displaces
`P2-G` (max/max, 15 ep)** or closes the classical gap. Cascade remains **not**
≥0.03 tissue F1 ahead of `C-mvalue-enet-G`.

### RBS-only frozen diagnostic (`rbs_linear_probe`)

**It did run** on all four screen cascade arms (inline; `rbs_enet` deferred with
screen `mbs_enet` policy). It was **missing from this report** until a bugfix:
`_metric_from_fold` omitted `rbs_*` from evaluation keys, so means/`task_comparison`
dropped those rows. **`P2-G` / `P4-G` still lack `rbs_linear_probe`** in metrics
(historical runs; only fold 0 has `all_gene_rbs.zarr` for P2).

| Arm | RBS linear F1 | MBS linear F1 | RBS age MAE | MBS age MAE | RBS sex AUROC | MBS sex AUROC |
|-----|-------------:|-------------:|------------:|------------:|--------------:|--------------:|
| `N-cascade-scalar-mean-max` | 0.385 | 0.375 | 14.67 | 14.38 | 0.773 | 0.711 |
| `N-cascade-scalar-max-mean` | 0.383 | 0.367 | 13.51 | 12.25 | 0.811 | 0.800 |
| `N-cascade-vector-mean-max` | 0.368 | 0.360 | **10.46** | **15.05** | **0.842** | **0.697** |
| `N-cascade-vector-max-max` | 0.329 | 0.367 | 12.10 | 14.86 | 0.803 | 0.687 |

**Read:** tissue is near-neutral across gene pooling. On **vector-mean-max**, RBS
beats MBS on age (~4.6 y) and sex (+0.15 AUROC) — evidence that **region→gene
pooling discards age/sex signal** even though vector e2e did not beat scalar P2.
Scalar arms show little age loss at the gene pool (sometimes MBS age is slightly
better). Classical age (8.15) remains far ahead of both RBS and MBS.

### Lock (unchanged)

Keep **`P2-G` max/max 15 ep**. No Tier-2 promotions. Next: Stage B GPU; optional
post-hoc `mbs_enet` / `rbs_enet` and backfill P2-G RBS probe if desired.
## Architecture questions (Stage A screen)

1. **CpG → region pool (mean vs max):** `mean-max` tissue F1=0.331 vs `max-max` 0.373; age MAE 20.713 vs 15.637. Prefer **`P2-G`** on this slice (check Pareto).
2. **Region → gene pool (mean vs max):** `max-mean` tissue F1=0.359 vs `max-max` 0.373; age MAE 21.356 vs 15.637. Prefer **`P2-G`** on this slice (check Pareto).
3. **Does scalar RBS discard information?** Vector arm `0.337` tissue vs P2 `0.373`; if vector does not beat scalar on age/sex, bottleneck is elsewhere.
4. **Gene pooling vs RBS:** `N-cascade-vector-mean-max` `rbs_*` tissue F1=0.368, age MAE=10.459, sex AUROC=0.842; same-arm MBS probe tissue=0.360, age=15.051, sex=0.697. Gene pooling is near-neutral on tissue; **age/sex often better on RBS** (pre–gene-pool), so some phenotype signal is lost at region→gene. Classical enet age MAE=8.150 remains the age ceiling.
5. **One-hop vs cascade:** One-hop `N-light-gene-max` tissue=0.122 / age=21.888 vs P2-G 0.373 / 15.637.
6. **One-scalar-per-gene bottleneck:** Gene aggregation still trails classical on age/sex; one scalar MBS/gene is **not yet adequate** unless a screen arm closes the gap.
7. **Best performance/compute:** Prefer landed P2/P4 (15 ep) over P5; promote Tier-1 (5 ep) arms only when Pareto/near-best, then confirm at 15 ep.

## Cascade arms (gene-linked CpGs only)

Primary **`mbs_e2e`** (test split only); **`mbs_linear_probe`** and **`mbs_enet`** are readouts of the **same frozen MBS**; **`rbs_*`** use gene-linked RBS. Contaminated pre-fix **`mbs_e2e`** shown as *invalid*.

| Arm | mbs_e2e F1 | linear probe F1 | mbs_enet F1 | age MAE (e2e) | sex AUROC (probe) | folds |
|-----|-----------:|----------------:|------------:|--------------:|------------------:|------:|
| P2-G | 0.373 (±0.038) | 0.373 | 0.385 (±0.053) | 15.637 | 0.784 | 3 |
| P4-G | 0.370 (±0.059) | 0.379 | 0.379 (±0.059) | 20.380 | 0.800 | 3 |
| N-cascade-scalar-max-mean | 0.359 (±0.075) | 0.367 | — | 21.356 | 0.800 | 3 |
| P5-G-max | 0.356 (±0.042) | 0.371 | — | 21.402 | 0.714 | 3 |
| N-cascade-vector-max-max | 0.343 (±0.063) | 0.367 | — | 21.458 | 0.687 | 3 |
| N-cascade-vector-mean-max | 0.337 (±0.036) | 0.360 | — | 22.753 | 0.697 | 3 |
| N-cascade-scalar-mean-max | 0.331 (±0.020) | 0.375 | — | 20.713 | 0.711 | 3 |
| N-light-gene-mean | 0.126 (±0.064) | 0.345 | 0.278 (±0.000) | 23.082 | 0.770 | 3 |
| N-light-gene-max | 0.122 (±0.102) | 0.295 | 0.131 (±0.000) | 21.888 | 0.731 | 3 |

## Classical arms (-G panel)

Same **51,375 gene-linked CpGs** as neural arms (`explicit_only`): ridge, elastic-net, HGB, PCA-SVA+ridge.

| Arm | tissue F1 | age MAE | age R² | sex AUROC |
|-----|----------:|--------:|-------:|----------:|
| C-mvalue-enet-G | 0.388 (±0.018) | 8.150 | 0.797 | 0.882 |
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

> Values populated only when `stage_a_per_epoch_eval: true` and `repr_diagnostics` logged on `mbs_e2e` (often empty for these short runs).

| Arm | Gene-score SD | Saturation frac | Const-score frac | Corr w/ mean-M |
|-----|:-------------:|:---------------:|:----------------:|:--------------:|
| A0 | — | — | — | — |
| A1 | — | — | — | — |
| A2 | — | — | — | — |
| A3 | — | — | — | — |
| A4/A7 | — | — | — | — |
| N0 | — | — | — | — |
| N1 | — | — | — | — |
| N2 | — | — | — | — |
| N3 | — | — | — | — |

## Parallel / follow-on work

- **Stage A required GPU arms** (`P2-G`, `P4-G`, `P5-G-max`, `C-mvalue-*-G`) are complete on `explicit_only`. Optional `P5-G-mean` was not run.
- **Stage A screen (sequential):** train one arm at a time and regenerate this report after each. Order: `N-light-gene-max` → `N-light-gene-mean` → mixed scalar cascades → vector cascades; promote Tier-2 (15 ep) only if Pareto/near-best.
- **Encoder parity (optional):** FlatDeepSet + HierarchicalDeepSet on same `gene_cols` if cascade does not lead classical by ≥0.03 F1.
- **Stage B (after lock):** fold-safe `C-mvalue-enetS`, `N-cascade-S`, `N-light-type`, `direct_cpg.zarr`, full-model fusion arms.

## Next

- **Stage A screen (sequential):** continue remaining Tier-1 arms after each landed light arm updates this report.
- Stage B (after lock): fold-safe `C-mvalue-enetS`, `N-cascade-S`, `N-light-type` (FlatDeepSetRegion), `N-mbs-posthoc-full-fusion` / `N-mbs-posthoc-mbs-direct`, plus `direct_cpg.zarr`.
- Milestone **7** 5×6 OOF remains blocked until Stage B completes.
