# Frozen P2-G probes — disease / cancer / BMI

Generated `2026-09-08T13:32Z` from
`scripts/run_7h_disease_cancer_frozen_probes.py` on
`stage0-7h-nine-pack-P2-G` (per-fold checkpoint, no encoder retrain).
Split: `hub-nine-pack-3fold-v1`. Raw: [`results.json`](results.json).

Logistic (binary) / ridge (BMI) on frozen MBS or RBS.

| Trait | n (usable) | MBS 3-fold | RBS 3-fold |
|-------|-----------:|------------|------------|
| **Cancer** case/control | 9 077 | AUROC **0.954** | AUROC **0.958** |
| **Disease** (any vs control) | 12 194 | AUROC **0.586** | AUROC **0.585** |
| **BMI** | 2 070 | MAE **9.87**, R² &lt; 0 | MAE **9.33**, R² ≈ 0 |

**Interpretation.** Cancer signal is strong on a frozen encoder that never
saw a cancer head. Broad disease is only modest — the pack label mixes 28
diagnoses. BMI is not yet a useful freeze-reuse target (negative R²).
Individual-disease probes: Alzheimer’s 3-fold AUROC **0.838** (see
`subtype_and_homogeneity.json`).
