# Stage 0 Milestone 7G — Methylation-only full evaluation

Status: **done** when this report exists with Ranking + winner.

## Budget and remaining ceilings

- Neural / classical loci: **65536** (matrix has **482379**).
- Epochs: **15**; restarts: **1**.
- Remaining ceiling: full **482 379** loci and a 2nd restart — deferred.
- Trees: sklearn **HistGradientBoosting** (not LightGBM).
- SVA: **10 train-fold PCA** surrogate variables (not Bioconductor `sva`).
- Metadata-only omitted from ranking (7E′ leakage alarm only).

## Topology under test

```text
CpG → typed region (gene | RBS) → RBS
  ├─ allocated to gene → MBS
  └─ orphan RBS
leftover CpG (no typed region / former TBS) → direct
late fusion: [orphan RBS | MBS | direct] → linear heads
```

No TBS arm (ADR 0009).

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
| `T-mean-gene` | Presence-aware **gene-mean** M-values → fold-fitted linear multitask heads. |
| `T-mean-region` | Presence-aware **region-mean** M-values → linear heads (strong tissue baseline). |
| `T-enet` | Transparent elastic-net on gene- or region-aggregated M-value features. |
| `C-mvalue-ridge` | All prefix loci → M-value → train-fold impute/scale → ridge (age) / SGD-L2 logistic (tissue, sex). |
| `C-mvalue-enet` | All prefix loci → M-value → elastic-net regression (age) / logistic enet (tissue, sex); **7G bake-off tissue winner** (F1≈0.334 on 65k). |
| `C-mvalue-hgb` | M-value → HistGradientBoosting (stands in for LightGBM) per trait. |
| `C-mvalue-sva` | M-value → PCA surrogate components (Bioconductor `sva` stand-in) → ridge/logistic. |
| `N-cascade-l1` | **7F/7G product topology:** CascadeDeepSet (CpG→RBS→MBS) + fold-fitted direct enet + **late fusion** on saved blocks; Level-1 z on direct branch. |

## Ranking (methylation-input only)

Tissue macro-F1 / balanced accuracy exclude classes with zero training examples in that fold (study-grouped folds can hold an entire rare tissue class out of train; no model can predict those by construction, so counting them would penalize every arm for a fold-construction artifact rather than model quality). `n classes scored` / `n excluded (zero-shot)` make the denominator explicit per arm.

| Arm | Tissue macro-F1 | Balanced acc | n classes scored | n excluded (zero-shot) | Sex AUROC | Age MAE | Age R² |
|-----|-----------------|--------------|-------------------|------------------------|-----------|---------|--------|
| T-mean-gene | 0.323 | 0.354 | — | 0 | 0.880 | 8.909 | 0.755 |
| T-mean-region | 0.330 | 0.359 | — | 0 | 0.908 | 7.946 | 0.802 |
| T-enet | 0.334 | 0.365 | — | 0 | — | 19.153 | 0.125 |
| C-mvalue-ridge | 0.288 | 0.313 | — | 0 | 0.910 | 6.372 | 0.859 |
| C-mvalue-enet | 0.334 | 0.377 | — | 0 | 0.894 | — | — |
| C-mvalue-hgb | 0.077 | 0.102 | — | 0 | 0.943 | 8.756 | 0.769 |
| C-mvalue-sva | 0.301 | 0.326 | — | 0 | 0.852 | 12.498 | 0.108 |

## Winner (Milestone 7 topology)

**`C-mvalue-enet`** — max mean tissue macro-F1, then min mean age MAE among methylation-input methods. Tissue macro-F1=0.334; sex AUROC=0.894; age MAE=—.

## ROC

Sex AUROC and tissue one-vs-rest curves in `figures/` come from **neural fusion** scores (`N-cascade-l1`), not only HGB.

## Artifacts

- `summary.json` — full dump including sex metrics
- `classical_baselines.json`
- `transparent_baselines.json`
- `arm_means.json`
- `figures/roc_tissue_ovr_fusion.png`, `figures/roc_sex_fusion.png`
