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

## Ranking (methylation-input only)

| Arm | Tissue macro-F1 | Balanced acc | Sex AUROC | Age MAE | Age R² |
|-----|-----------------|--------------|-----------|---------|--------|
| N-cascade-l1 | 0.093 | 0.125 | 0.911 | 8.441 | 0.782 |
| T-mean-gene | 0.323 | 0.354 | 0.880 | 8.909 | 0.755 |
| T-mean-region | 0.330 | 0.359 | 0.908 | 7.946 | 0.802 |
| T-enet | 0.334 | 0.365 | — | 19.153 | 0.125 |
| C-mvalue-ridge | 0.288 | 0.313 | 0.910 | 6.372 | 0.859 |
| C-mvalue-enet | 0.334 | 0.377 | 0.894 | — | — |
| C-mvalue-hgb | 0.077 | 0.102 | 0.943 | 8.756 | 0.769 |
| C-mvalue-sva | 0.301 | 0.326 | 0.852 | 12.498 | 0.108 |

## Winner (Milestone 7 topology)

**`C-mvalue-enet`** — max mean tissue macro-F1, then min mean age MAE among methylation-input methods. Tissue macro-F1=0.334; sex AUROC=0.894; age MAE=—.

## ROC

Sex AUROC and tissue one-vs-rest curves in `figures/` come from **neural fusion** scores (`N-cascade-l1`), not only HGB.

## Follow-up

Cascade tissue macro-F1 (~0.09) vs transparent/classical (~0.33): see
[`docs/plans/milestone-7g-cascade-tissue-investigation.md`](../../../docs/plans/milestone-7g-cascade-tissue-investigation.md).

## Artifacts

- `summary.json` — full dump including sex metrics
- `classical_baselines.json`
- `transparent_baselines.json` (local only, gitignored; regenerate from driver)
- `arm_means.json`
- `figures/roc_tissue_ovr_fusion.png`, `figures/roc_sex_fusion.png`
