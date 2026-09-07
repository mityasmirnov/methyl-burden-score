# Plan: 7G″ expression-auxiliary pilot (deferred)

> **Canonical milestone: 13** (deferred; **after** Milestone **12** OOF).
> Historical **7G″**. See [`MILESTONE_INDEX.md`](MILESTONE_INDEX.md) and
> [`milestone-13-expression-auxiliary.md`](milestone-13-expression-auxiliary.md).


Status: **deferred** — not a gate for Milestones **10** / **11** / **12** OOF.
Runs **after** final OOF (or an explicit interim pretrained checkpoint).
Parent: [`milestone-7g-prime-matched-probe-lightweight.md`](milestone-7g-prime-matched-probe-lightweight.md).

## Scope

Test whether **gene-expression auxiliary supervision** (continue-training or
finetuning the methylation encoder) teaches more biologically interpretable
promoter/body aggregation than phenotype supervision alone. Inspired by
RSMethy-Net (CSBJ 0138), treated as a **hypothesis generator only** — not an
architecture to copy.

## Data download (required before GPU)

1. Download matched expression matrices (RNA-seq and/or expression arrays) for
   cohorts that overlap canonical Hub / TCGA methylation samples.
2. Store under `$MBS_DATA_ROOT` with manifests, checksums, and study/sample IDs.
3. Publish an overlap census (n samples, n genes, platforms) under
   `reports/inspection/` before any expression train/finetune job.

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

1. Complete Milestones **9–11** (gene-only selection, scale, fold panel).
2. Complete Milestone **12** final OOF (or document an interim **10d**
   checkpoint if OOF is postponed by ADR).
3. **Download + catalog expression** (see above) and overlap census.
4. Run 7G″ / Milestone **13** continue-train and/or finetune on a **bounded**
   RNA-overlap subset (not full Hub OOF).
5. If expression amends the encoder for product use, document the new
   checkpoint contract — do **not** silently replace OOF scores.

## Non-goals

- Replace CascadeDeepSet with RSMethy-Net per-gene encoders.
- Block **10** / **11** / **12** on 7G″ completion.
- Start expression GPU before download + overlap census.

## Done when (future)

- Expression download + manifests on disk; overlap census reported.
- Report under `reports/inspection/stage0_7g_double_prime_expression/` with
  continue-train / finetune arms, residualization ablation, and phenotype +
  expression metrics on study-held-out folds.
