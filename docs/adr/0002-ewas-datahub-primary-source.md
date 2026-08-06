# ADR 0002: EWAS Data Hub as primary Stage 0 open data source

## Status

Accepted

## Context

Stage 0 originally inspected a tiny CpGCorpus GSE/GPL subset
(`GSE125367` / `GPL21145`) and documented a broader labeling cohort against
requester-pays S3. Many of those labeling GSEs are absent from the public
CpGCorpus prefix, while all are present as per-sample beta files under CNCB
EWAS Data Hub `EWAS_db/`. Data Hub also publishes large GMQN-normalized baseline
packs (tissue, blood, brain, disease, age, …) under `download/`.

Relying on paywalled or incomplete corpora restricts reproducibility and
complicates pilot matrix construction. Download and layout for Atlas + Data Hub
are already documented in [`docs/EWAS_DATA.md`](../EWAS_DATA.md).

## Decision

1. **Primary open source** for Stage 0 pilot canonical matrices and subsequent
   open-scale training is the CNCB **EWAS Data Hub** (`EWAS_db/` and/or small
   baseline subsets under `download/`).
2. **EWAS Atlas** remains the curated association resource for later enrichment
   / knowledge-graph validation (post–Stage 0).
3. **CpGCorpus** is an optional alternate path (requester-pays). Milestone 1
   evidence based on CpGCorpus inspection stays valid and is not re-opened.
4. **Model milestone order is unchanged:** annotation graph → static locus
   features → pilot matrix → flat DeepRVAT baseline → hierarchical model →
   study-grouped cross-fitting. PROTRIDER-style autoencoders, ComBat-met, and
   REGENIE export are deferred until after that core path
   ([`docs/STRATEGIC_PLAN.md`](../STRATEGIC_PLAN.md)).

## Consequences

- [`docs/TODO_PIPELINE.md`](../TODO_PIPELINE.md) milestone 4 and experiment
  protocol pilot wording target EWAS Data Hub, not CpGCorpus Arrow, as the
  default.
- [`docs/CPGCORPUS_STAGE0.md`](../CPGCORPUS_STAGE0.md) documents an optional /
  historical path; agents must not treat CpGCorpus sync as required for
  milestones 2–7.
- Existing `make download-ewas-datahub` / Atlas scripts remain the supported
  acquisition path; greenfield FTP-only scaffolds are out of scope.
- No change to Deep Set architecture contracts in
  [`docs/ARCHITECTURE.md`](../ARCHITECTURE.md).
