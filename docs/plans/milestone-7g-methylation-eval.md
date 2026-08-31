# Plan: Milestone 7G — Methylation-only full evaluation

Status: **done** (closed 2026-08-27). Follow-up:
[`milestone-7g-cascade-tissue-investigation.md`](milestone-7g-cascade-tissue-investigation.md).
Report target: [`reports/inspection/stage0_7g_methylation_eval/`](../../reports/inspection/stage0_7g_methylation_eval/).
Normative: [ADR 0006](../adr/0006-multipath-noncoding-scores.md),
[ADR 0007](../adr/0007-crossfit-prerequisites.md),
[ADR 0008](../adr/0008-score-identifiability.md),
[ADR 0009](../adr/0009-drop-tbs-scores.md),
[TODO_PIPELINE.md](../TODO_PIPELINE.md) §7G,
[post-v0-scientific-programme.md](post-v0-scientific-programme.md) §7G.
Depends on: 7F ([`milestone-7f-rbs-gene-direct.md`](milestone-7f-rbs-gene-direct.md)).

## Scope and acceptance

**Done when** (from TODO):

- 7F cascade retrained on frozen `hub-ats-7e-3fold-v1` with budget above
  7E’s 2-epoch / 8 192-locus ceiling (or an explicit remaining ceiling)
- Classical M-value models use methylation matrices only
- ROC for **sex** and **tissue one-vs-rest** from neural fusion scores
  (age stays MAE / R², no ROC)
- Sex metrics present in the summary dump; region-mean transparent arm named
  if kept (`T-mean-region`)
- Report under `reports/inspection/stage0_7g_methylation_eval/` names the
  Milestone 7 topology from **methylation-input methods only**

Do **not** overwrite the 7F report, v0.1 freezes, or `deepmat-data-age-tissue-sex-v1`.

## Locked decisions

| Choice | Decision | Why |
|--------|----------|-----|
| Splits | Frozen `hub-ats-7e-3fold-v1` only | Same studies/folds as 7E |
| Neural budget | **65 536 loci, 15 epochs, 3 folds × 1 restart** | 8× loci and 7.5× epochs vs 7E |
| Remaining ceiling | 482 379 loci; no 2nd restart | Named honestly in the report |
| Classical loci | Same **65 536** prefix as neural | Matched comparison |
| Trees | sklearn HistGradientBoosting (no LightGBM) | Same family; no new dep |
| PCA-SVA | 10 train-fold PCs (not Bioconductor `sva`) | Close 7E gap without rpy2 |
| Encoder | 64 / GELU / dropout 0.1 / LayerNorm | Parameter-matched |
| Fusion | Saved `[orphan RBS \| MBS \| direct]` → linear heads | Not region-mean tables |
| Metadata-only | Omit from ranking | Leakage alarm from 7E′ only |
| Task masks | Honor age/tissue/sex masks | Unlabeled ≠ class 0 |
| Disconnect | `nohup` + per-fold skip-if-done | Survive session disconnect |

## Schemas / contracts

Config: `configs/experiment/stage0_7g_methylation_eval.yaml`.

Score export under `$MBS_ARTIFACT_ROOT/runs/<run_id>/fold_<i>/scores/`
(same layout as 7F; no `tbs.zarr`).

Report:

```text
reports/inspection/stage0_7g_methylation_eval/
├── analysis.md
├── summary.json
├── classical_baselines.json
├── arm_means.json
└── figures/
```

Named cells (methylation input only):

| Arm | Family |
|-----|--------|
| `N-cascade-l1` | neural (7F topology, saved-score fusion) |
| `T-mean-gene` | transparent |
| `T-mean-region` | transparent |
| `T-enet` | transparent |
| `C-mvalue-ridge` | classical |
| `C-mvalue-enet` | classical |
| `C-mvalue-hgb` | classical |
| `C-mvalue-sva` | classical |

## Data / artifact flow

```text
hub-ats-7e-3fold-v1 + matrix-hub-age-tissue-sex-full-v1 + graph-v2
  → CascadeDeepSet (65536 loci, 15 epochs) + Level-1 direct enet
  → write score Zarrs (orphan RBS | MBS | direct)
  → late fusion + sex AUROC / tissue OvR from fusion proba
  → classical M-value ridge / enet / HGB / PCA-SVA (same prefix)
  → T-mean-gene / T-mean-region / T-enet
  → reports/inspection/stage0_7g_methylation_eval/
```

CLI:

```bash
source scripts/activate_data_environment.sh
bash scripts/train_7g_methylation_eval_background.sh   # nohup; survives disconnect
bash scripts/status_7g_methylation_eval.sh
```

## Non-goals / deferred

- Full 482 379-locus neural or classical
- LightGBM package; Bioconductor `sva`; 3×2 restarts
- Metadata-only in ranking; TBS scores; Milestone **7** 5×6 OOF
- Overwriting 7F report or v0.1 freezes

## Outcome (2026-08-27)

| Arm | Tissue macro-F1 | Sex AUROC | Age MAE | Age R² |
|-----|-----------------|-----------|---------|--------|
| **C-mvalue-enet** (ranking winner) | **0.334** | 0.894 | — | — |
| T-mean-region | 0.330 | 0.908 | 7.95 | 0.802 |
| N-cascade-l1 (7F topology) | 0.093 | 0.911 | 8.44 | 0.782 |

Ranking rule (max tissue F1, then min age MAE) selects a **classical** arm, not
the neural cascade. Milestone **7** OOF may still export 7F product scores; see
tissue probe plan before claiming cascade parity on tissue.

## Open questions

Resolved for 7G scope. Tissue-head gap →
[`milestone-7g-cascade-tissue-investigation.md`](milestone-7g-cascade-tissue-investigation.md).
Milestone **7** topology vs ranking winner may need an ADR.
