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
| `P4-G` | Like P2-G but **mean/mean** pooling; same gene-only panel and **`mbs_e2e`** metric. |
| `P2-G` | CascadeDeepSet on **gene-linked CpGs only** (`gene_cols`); max/max pooling; 15 epochs; P2 loss weights; primary metric **`mbs_e2e`** (end-to-end MBS heads). |
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
| `N-full` | **Stage B full model:** late fusion **`fusion_full`** = orphan RBS + MBS + direct columns. |
| `N-mbs-direct-only` | **Orphan ablation:** late fusion **`fusion_mbs_direct`** = MBS + direct only (no orphan RBS block). |

### Evaluation modes (cascade metrics)

- **`mbs_e2e`**: End-to-end **`MultitaskHeads`** on MBS only (no orphan/direct columns) — **Stage A primary metric**.
- **`mbs_linear_probe`**: CPU linear probe fit on **saved MBS matrix** only (representation check, not fusion).
- **`fusion_full`**: Late fusion on **`[orphan_rbs | mbs | direct_contrib]`** columns.
- **`fusion_mbs_direct`**: Late fusion on **`[mbs | direct_contrib]`** — orphan RBS ablation.

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

## Cascade arms (gene-linked CpGs only)

Primary **`mbs_e2e`**; secondary **`mbs_linear_probe`**. Age/sex from **`mbs_e2e`** heads.

| Arm | mbs_e2e F1 | linear probe F1 | age MAE | age R² | sex AUROC/F1 | folds |
|-----|-----------:|----------------:|--------:|-------:|-------------:|------:|
| P4-G | 0.695 (±0.054) | 0.374 | 11.961 | 0.554 | 0.778 | 3 |
| P2-G | 0.682 (±0.030) | 0.380 | 10.128 | 0.667 | 0.807 | 3 |
| P5-G-mean | 0.672 (±0.029) | 0.368 | 13.438 | 0.456 | 0.810 | 3 |
| P5-G-max | 0.669 (±0.052) | 0.366 | 11.370 | 0.588 | 0.790 | 3 |

## Classical arms (-G panel)

Same **56 214 gene-linked CpGs** as neural arms: ridge, elastic-net, HGB, PCA-SVA+ridge.

| Arm | tissue F1 | age MAE | age R² | sex AUROC |
|-----|----------:|--------:|-------:|----------:|

## Locked architecture (Stage B input)

- **Cascade arm:** `P4-G`
- **Pooling (CpG / region):** `mean` / `mean`
- **Epoch ceiling:** 15
- **Best classical (-G):** `None`
- **Cascade clearly ahead (≥0.03 tissue F1):** None

Cascade leads classical `-G` by ≥0.03 F1; optional Flat/Hier parity runs are **not** required.

Orphan ablation arm **P2-orphan-ablation** not run.

## Parallel / follow-on work

- **Stage A finish:** classical `-G` (CPU) + orphan ablation (GPU) — can run in parallel.
- **Encoder parity (optional):** FlatDeepSet + HierarchicalDeepSet on same `gene_cols` if cascade does not lead classical by ≥0.03 F1.
- **Stage B (after lock):** fold-safe `C-mvalue-enetS`, `N-cascade-S`, `N-light-type`, `direct_cpg.zarr`, full-model fusion arms.

## Next

- Stage B: fold-safe `C-mvalue-enetS`, `N-cascade-S`, `N-light-type` (FlatDeepSetRegion), `N-full` / `N-mbs-direct-only`, plus `direct_cpg.zarr`.
- Milestone **7** 5×6 OOF remains blocked until Stage B completes.
