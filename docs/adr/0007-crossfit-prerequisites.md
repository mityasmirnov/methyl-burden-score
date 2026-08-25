# ADR 0007: Insert Milestones 7A–7E before final OOF cross-fitting

## Status

Accepted

## Context

Milestones 1–6 are done. Flat deepMAT-v0.1 and hierarchical-v0.1 on the
age/tissue/sex GSM union are useful baselines; hierarchical underperformed
flat on tissue accuracy and age MAE. The repository previously treated
Milestone 7 (5-fold × up to 6-restart study-grouped OOF scoring) as the
immediate next gate.

Launching that expensive protocol on the current architecture would lock in:

- fragmented catalog / incomplete nine-pack matrices and multi-label disease;
- sample-count-only splits;
- one-scalar residual noncoding path;
- inconsistent head centering and missing seed/task-balance machinery;
- no fold-fitted robust normalization channel.

## Decision

1. Amend Stage 0 order after Milestone 6
   ([ADR 0002](0002-ewas-datahub-primary-source.md) /
   [ADR 0003](0003-milestone-5b-phenotype-registry.md)):

   ```text
   … → hierarchical (6) →
   7A harmonized release + phenotype census →
   7B complete nine-pack matrices →
   7C architecture corrections (heads, splits, multi-path scores) →
   7D fold-fitted normalization (Level 1 required) →
   7E development CV (3 folds × 2 restarts) →
   7 final OOF cross-fitting (5 folds × up to 6 restarts) →
   8 optional layers
   ```

2. Milestone **7** remains the OOF score-matrix deliverable. Its status stays
   `pending` and is **blocked until 7A–7E** acceptance criteria are met.
3. A small 3-fold / 1-restart smoke of *existing* train machinery is allowed
   for plumbing; it does **not** complete Milestone 7, must not overwrite
   frozen v0.1 runs, and is **not** development CV (7E).
4. Freeze (do not overwrite) reference artifacts:
   - `deepMAT-flat-v0.1` ← `stage0-flat-deeprvat-age-tissue-sex-full-v1`
   - `deepMAT-hierarchical-v0.1` ← `stage0-hier-deeprvat-age-tissue-sex-full-v1`
   - `deepmat-data-age-tissue-sex-v1` ← `matrix-hub-age-tissue-sex-full-v1` (+ phenotype table)

## Consequences

- [`TODO_PIPELINE.md`](../TODO_PIPELINE.md) orders **7A → 7B → 7C → 7D → 7E → 7**.
  After 7A landed the gate was **7B**. As of 2026-08-25, **7A–7D are done**;
  the current gate is **7E** (graph-v2 on disk; multi-path unblocked for
  graph-v2).
- Build brief:
  [`plans/post-v0-scientific-programme.md`](../plans/post-v0-scientific-programme.md).
- Storage / catalog: [ADR 0005](0005-catalog-matrix-independence.md).
- Noncoding paths: [ADR 0006](0006-multipath-noncoding-scores.md).
- Score orientation / predictive vs constraint: [ADR 0008](0008-score-identifiability.md).
- EWAS_db All-Data mirror progress is **not** a gate for 7A–7E.
- **Do not retrain v0.1.** After 7A the first coding deliverable was **7B**.
  7E (gene-only) starts only after 7B–7D, which are now closed. Milestone **7**
  stays blocked until 7E.

## Non-goals

- Reopening Milestones 1–6.
- Selecting architecture on reconstruction loss alone (7D).
- Claiming hierarchical-v0.1 as the preferred phenotype model.
