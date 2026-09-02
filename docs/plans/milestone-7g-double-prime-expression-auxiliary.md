# Plan: 7G″ expression-auxiliary pilot (deferred)

Status: **deferred** — not a gate for 7G′ Stage A/B or Milestone **7** OOF.
Parent: [`milestone-7g-prime-matched-probe-lightweight.md`](milestone-7g-prime-matched-probe-lightweight.md).

## Scope

Test whether **gene-expression auxiliary supervision** teaches more biologically
interpretable promoter/body aggregation than phenotype supervision alone. Inspired
by RSMethy-Net (CSBJ 0138), treated as a **hypothesis generator only** — not an
architecture to copy.

## Why not copy RSMethy-Net directly

- Random within-cancer CV; no cross-study generalization.
- TCGA / HM450 only; ±10 Mb gene windows.
- Per-gene fixed-size encoders — incompatible with probe-agnostic DeepSet MBS.
- Reported region-structure gain ~**0.02 median expression R²** in 5/6 cancers.

## Auxiliary loss (sketch)

For sample \(s\) and gene \(g\) with expression mask \(m_{s,g}\):

\[
\mathcal L =
\mathcal L_{\text{phenotypes}}
+
\lambda_{\text{expr}}
\sum_{s,g} m_{s,g}
\left[
\operatorname{SmoothL1}(\hat e_{s,g}, e_{s,g})
+
\gamma(1-\operatorname{PCC})
\right]
\]

Never impute missing expression targets; mask unknown RNA.

## Pilot arms

| Arm | Representation | Purpose |
|-----|----------------|---------|
| A | Locked phenotype P*-G (Stage A winner) | Reference |
| B | Expression-only shared encoder | Does expression teach aggregation? |
| C | Phenotypes + expression auxiliary | Primary candidate |
| D | Shared region-type embedding | Lightweight baseline |
| E | Shared trunk + promoter/body adapters | Paper-inspired, parameter-efficient |
| F | Signed promoter/body outputs | Test scalar MBS direction loss |

## Safeguards

- Study- or cancer-grouped splits (not random samples).
- Fold-fitted methylation normalization (Level-1 inside outer train).
- Raw vs tissue/study-**residualized** expression (avoid tissue-expression shortcut).
- Shared encoder across genes; no CpG/gene ID embeddings if agnosticism required.
- Platform downsampling: full → HM450 → EPIC overlap.
- Evaluate expression R²/PCC **and** phenotype metrics; ADR 0008 orientation anchor.
- Expression–methylation graph edges only from independent reference or outer-train fit.

## Sequencing

1. Complete honest **7G′ Stage A** lock (test-only `mbs_e2e`, `explicit_only`, matched `-G`).
2. Complete **7G′ Stage B** on fold-selected panels.
3. Run 7G″ pilot on a **bounded** RNA-overlap subset (not full Hub OOF).
4. If arm C amends the encoder, update Stage A lock **before** Milestone **7** OOF.

## Non-goals

- Replace CascadeDeepSet with RSMethy-Net per-gene encoders.
- Block 7G′ or Milestone **7** on 7G″ completion.
- Ingest TCGA RNA in this scaffolding change set.

## Done when (future)

- Report under `reports/inspection/stage0_7g_double_prime_expression/` with arms A–F,
  residualization ablation, and phenotype + expression metrics on study-held-out folds.
