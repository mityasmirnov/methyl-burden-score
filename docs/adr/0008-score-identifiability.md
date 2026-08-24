# ADR 0008: MBS orientation identifiability; predictive vs constraint scores

## Status

Accepted

## Context

deepMAT maps samples to gene-level (and, from 7C, region/tile) scores that feed
linear phenotype heads. With centered sigmoid MBS in `[0, 1]` and unconstrained
linear weights, the map

```text
MBS → 1 − MBS,    head weights → −weights
```

yields identical predictions. Fold- and restart-specific scores can therefore
have **incompatible orientation**. Averaging OOF MBS without an orientation
anchor is not scientifically defined.

Separately, a gene-level **predictive representation** (what deepMAT trains
today) is not a LOEUF-like **constraint** score. Constraint would be estimated
later from the tissue-stratified distribution and depletion of OOF residual
burdens in an appropriate reference population.

v0.1 residual-only ablations also evaluated the **first 512 ordered** holdout
samples rather than a stratified subset; that sampling choice is not evidence
about noncoding biology ([ADR 0006](0006-multipath-noncoding-scores.md)).

## Decision

1. Before averaging OOF scores across folds or restarts (Milestone **7**), define
   an **orientation anchor**. Preferred Stage 0 options:
   - separate **hyper-** and **hypomethylation** burden channels; or
   - a single burden **magnitude** anchored to absolute fold-fitted robust
     M-deviation (sign/orientation recorded in the score manifest).
2. Persist orientation metadata (`score_polarity`, anchor recipe, fold/restart
   IDs) in `score_manifest.json`.
3. Keep **predictive MBS/RBS/TBS** distinct from a future **methylation
   constraint** score. Do not advertise OOF MBS as constraint or LOEUF analogue.
4. Residual-path and branch ablations must use stratified (or full) holdout
   evaluation, not an ordered prefix of samples.

## Consequences

- Milestone **7C** implements the anchor in the scoring contract; Milestone **7**
  may not ensemble unaligned scores.
- [`ARCHITECTURE.md`](../ARCHITECTURE.md),
  [`DATA_CONTRACT.md`](../DATA_CONTRACT.md), and
  [`plans/post-v0-scientific-programme.md`](../plans/post-v0-scientific-programme.md)
  record the contract.
- Constraint-score estimation is **§8 / post–Stage 0**.

## Non-goals

- Estimating constraint scores in 7A–7E.
- Changing frozen v0.1 checkpoints.
