# Plan: 7G cascade tissue-head investigation (post-bake-off)

Status: **in progress** (P0–P3 complete 2026-08-31; P4/P5 implementation
landed, three-fold GPU results pending).
Parent: [`milestone-7g-methylation-eval.md`](milestone-7g-methylation-eval.md),
[`reports/inspection/stage0_7g_methylation_eval/analysis.md`](../../reports/inspection/stage0_7g_methylation_eval/analysis.md).
Normative: [ADR 0007](../adr/0007-crossfit-prerequisites.md),
[ADR 0008](../adr/0008-score-identifiability.md),
[ADR 0009](../adr/0009-drop-tbs-scores.md).

## Problem statement

7G closed the methylation-only bake-off on frozen `hub-ats-7e-3fold-v1` at
**65 536 loci / 15 epochs / 3×1**. The **7F product topology** (`N-cascade-l1`)
achieved strong **sex** (AUROC 0.911) and reasonable **age** (MAE 8.4 y, R²
0.78) but **weak tissue** (macro-F1 **0.093**). Transparent region means and
classical M-value elastic-net reached **~0.33–0.34** tissue macro-F1 on the
same folds.

Before launching Milestone **7** OOF with the cascade as the shipped encoder,
we need a bounded investigation: is the gap fixable by hyperparameters and
tissue-head design, or is gene-level aggregation intrinsically lossy for tissue
at this budget?

**This plan does not block Milestone 7** if programme leadership accepts
classical M-value enet as the methylation comparator winner and still runs OOF
for **product scores** (MBS / orphan RBS / direct) on the 7F topology. It
**does** block claiming the cascade “won” tissue without further evidence.

## Scope and acceptance

**Done when:**

- A short inspection report under
  `reports/inspection/stage0_7g_cascade_tissue_probe/` documents:
  - At least **three** ablation arms (see table below) on the **same** frozen
    folds, with tissue macro-F1, sex AUROC, and age MAE/R² per fold and means.
  - A **diagnosis section**: which failure mode best explains 7G (fusion
    bottleneck vs end-to-end head vs aggregation vs loss balance vs capacity).
  - A **recommendation** for Milestone 7: proceed cascade OOF as-is, run a
    narrowed hyperparameter grid first, or record an ADR that product OOF uses
    7F topology for **scores** while tissue prediction may remain a separate
    classical arm for benchmarking.
- Unit or integration smoke: one fixture fold proves each new code path runs.
- No overwrite of `stage0_7g_methylation_eval/` or v0.1 freezes.

**Out of scope for this probe:** full 482 379 loci, 5×6 OOF, LightGBM, ComBat,
architecture changes that drop RBS→gene/direct (ADR 0009).

## Locked decisions

| Choice | Decision | Why |
|--------|----------|-----|
| Splits | Frozen `hub-ats-7e-3fold-v1` only | Match 7G comparators |
| Base budget | **65 536 loci, 15 epochs, 1 restart** | Same as 7G unless an arm explicitly tests “more epochs” |
| Encoder | 64 / 32 / GELU / dropout 0.1 / LayerNorm | Hold architecture constant; vary heads and losses |
| Product topology | RBS→gene + direct; no TBS | ADR 0009 |
| Metrics | Tissue macro-F1 (primary), sex AUROC, age MAE/R² | Same as 7G ranking |
| Task masks | Honor age/tissue/sex masks | Scientific invariant |
| Selection | Report val-checkpoint tissue macro-F1 per fold | Align training with tissue goal |

## Ablation arms (minimum set)

Run on all three outer folds; column names are report labels.

| Arm | What changes | Hypothesis |
|-----|--------------|------------|
| **P0-baseline** | Replay 7G `N-cascade-l1` config (sanity) | Reproduce ~0.09 tissue F1 |
| **P1-fusion-tissue-heavy** | Late fusion only: `LogisticRegression`/`SGD` with `class_weight=balanced` on saved scores; optional PCA(32) on score block | Linear fusion underfits multi-class tissue |
| **P2-end2end-tissue-weight** | Cascade train: `tissue_loss_weight=3.0`, `age_loss_weight=0.3`, `sex_loss_weight=1.0`; checkpoint on val tissue macro-F1 | Age/sex gradients dominate equal weights |
| **P3-region-head-bypass** | Add eval arm: tissue head on **region-mean** features built from the same 65k prefix (no extra train) | Gene max-pool destroys tissue signal |
| **P4-pooling-mean** | Cascade: region/gene pooling **mean** instead of max (config flag) | Max pool keeps age/sex drivers, drops tissue diversity |
| **P5-epochs-30** | P2 settings but **30 epochs** + early stop on val tissue F1 | 15 epochs insufficient for tissue head |

Optional if P2–P5 fail to reach **0.20** mean tissue F1: **P6-direct-heavy-fusion**
(drop MBS block from fusion input) to test whether gene aggregation is the
bottleneck.

## Hyperparameter grid (narrow)

Only after P2/P5 identify sensitivity; do not full-grid search.

| Parameter | 7G value | Probe values |
|-----------|----------|--------------|
| `learning_rate` | 0.001 | 3e-4, 1e-3, 3e-3 |
| `tissue_loss_weight` | 1.0 | 0.5, 1.0, 3.0, 5.0 |
| `age_loss_weight` | 1.0 | 0.3, 1.0 |
| `dropout` | 0.1 | 0.0, 0.2 |
| `max_epochs` | 15 | 15, 30 (with early stop) |
| Fusion tissue solver | linear multitask | balanced logistic / SGD one-vs-rest |

Record best **per fold** on val tissue F1; report test metrics from that
checkpoint only (no test peeking).

## Schemas / contracts

Config sketch: `configs/experiment/stage0_7g_cascade_tissue_probe.yaml`
(fork of `stage0_7g_methylation_eval.yaml` with `arms:` list).

Report layout:

```text
reports/inspection/stage0_7g_cascade_tissue_probe/
├── analysis.md          # diagnosis + Milestone 7 recommendation
├── arm_means.json
├── per_arm/
│   ├── P0-baseline.json
│   └── ...
└── figures/
    └── tissue_f1_bars.png
```

Code touchpoints (expected minimal diffs):

- `src/mbs/training/cascade_loop.py` — loss weights, checkpoint criterion,
  pooling override, epoch early-stop hook.
- `src/mbs/training/late_fusion.py` — optional balanced multiclass solver for
  tissue block only.
- `scripts/run_7g_cascade_tissue_probe.sh` — idempotent runner (reuse 7G
  skip-if-done for P0 when scores exist).
- `scripts/write_7g_cascade_tissue_probe_report.py` — aggregate arms.

## Data / artifact flow

```mermaid
flowchart LR
  splits[hub-ats-7e-3fold-v1]
  mat[matrix-hub-age-tissue-sex-full-v1 prefix 65536]
  graph[graph-grch38-gencode38-cgi-tile-v2]
  splits --> train
  mat --> train
  graph --> assign[cascade_assign]
  assign --> train[CascadeDeepSet train]
  train --> scores[score Zarrs MBS RBS direct]
  scores --> fusion[late fusion ablations]
  train --> e2e[end-to-end multitask metrics]
  fusion --> report[stage0_7g_cascade_tissue_probe]
  e2e --> report
```

Reuse 7G fold score Zarrs for **P0** and **P1** when `skip-if-done` and
manifest hashes match.

## Non-goals / deferred

- Replacing Hub GMQN or Level-1 MAD (7D contract unchanged).
- Metadata-only arms; TBS scores; Milestone **7** full OOF in this probe.
- Claiming tissue parity with `C-mvalue-enet` without meeting the same masks
  and folds.

## Resolved findings and remaining checks

1. **Product vs phenotype:** the same phenotype-trained model exports MBS (and
   later orphan/direct features) **and** supports phenotype heads — not separate
   topologies. `C-mvalue-enet` (0.334 tissue F1) used all 65 536 prefix columns.
2. **P2 ~0.376 is not MBS-only evidence:** training supervises MBS; test metrics
   late-fuse `[orphan_rbs | mbs | direct_contrib]`. Corrected arms are **P2-G**
   … **P5-G** with **C-mvalue-enet-G** on identical gene-linked CpGs only.
3. **Checkpoint policy:** P2 best epochs 15/12/9 on validation tissue F1.
   Simulate patience 5 on stored val histories before locking; do not assume
   earlier stop is better (val F1 poorly tracked external-test F1).
4. **Study confounding:** per-fold composition tables remain interpretive only.

## Sequencing

1. P0–P3 remain **historical evidence** (committed report).
2. Implement gene-col filter + MBS-only evaluation mode.
3. Run **7G′ Stage A** (`P2-G`, `P4-G`, `P5-G`, `C-mvalue-enet-G`).
4. Run **7G′ Stage B** (fold-selected panel, full model, `direct_cpg.zarr`).
5. Start Milestone **7** 5×6 OOF.

## Phase 2 → superseded by 7G′ Stage A

The original P4/P5 grid (mean pooling, 30-epoch ceiling, fusion refit on saved
scores) measured **late-fusion** performance on the full 65 536 prefix. That
does not answer whether **gene-level aggregation** beats elastic-net on the
same gene-linked CpGs.

**Corrected immediate benchmark:** see
[`milestone-7g-prime-matched-probe-lightweight.md`](milestone-7g-prime-matched-probe-lightweight.md)
Stage A. P0–P3 and any uncorrected P4/P5 numbers are provisional until gene-only
arms complete.

During architecture selection, treat orphan-region CpGs as **non-gene excluded
input**. In the full model, qualified orphan regions stay one column per
`region_id`; unqualified regions route to direct.
