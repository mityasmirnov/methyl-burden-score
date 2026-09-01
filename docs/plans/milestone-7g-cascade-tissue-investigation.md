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

1. **Resolved by ADR 0010:** product score export remains deepMAT
   MBS/qualified orphan RBS/direct CpGs; phenotype comparators are ranked
   separately. `C-mvalue-enet` remains the current tissue comparator.
2. **Checkpoint policy resolved for P2:** selection is validation tissue
   macro-F1 then age MAE; best epochs were 15, 12 and 9. P5 adds an auditable
   early-stop record.
3. **Study confounding remains interpretive:** the report includes per-fold
   study/tissue composition; final claims remain held-out-study claims, not
   population prevalence claims.

## Sequencing

1. Close 7G in `TODO_PIPELINE.md` and commit inspection artifacts.
2. Implement P0–P3 (smallest code diff) + report writer.
3. Run probe on CPU/GPU via background script; idempotent folds.
4. Run P4/P5 and refusion grid; lock the Phase-2 winner.
5. Complete 7H before starting 5×6 OOF.

## Phase 2 (in progress)

P0–P3 completed; **mean tissue F1 0.376 (P2)** clears the 0.20 gate.
Implementation now includes:

- **P4**: P2 loss weights with mean CpG→region and region→gene pooling;
- **P5**: P2 loss weights, 30-epoch ceiling, validation-tissue macro-F1
  checkpointing, patience 8 and minimum improvement 0.001;
- narrow refusion cells: P2/P4/P5 saved scores with standard versus
  class-balanced logistic tissue fusion;
- strict three-fold completeness checks before an arm enters the report.

Performance remains **unknown** until all three P4 and P5 fold artifacts exist
on the data/GPU host. P0–P3 results must not be relabelled as P4/P5 evidence.

ADR 0010 separates product score export from phenotype comparators. The next
gate after Phase 2 is the fold-safe `C-mvalue-enetS`/deepMAT panel benchmark
and direct-CpG association export in
[`milestone-7h-fold-safe-probe-panel-benchmark.md`](milestone-7h-fold-safe-probe-panel-benchmark.md).
