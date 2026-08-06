# ADR 0003: Insert Milestone 5b (phenotype registry and multi-pack eval)

## Status

Accepted

## Context

ADR 0002 fixed the Stage 0 model path as: annotation graph → static features →
pilot matrix → flat DeepRVAT baseline → hierarchical model → study-grouped
cross-fitting. Milestone 5 (flat baseline) is complete on a single pilot study
(`GSE35069`). Starting the hierarchical model on that smoke pilot alone would
not show whether the shared scorer learns methylation biology versus
study/platform structure.

EWAS Data Hub already publishes phenotype-family baseline packs (age, tissue,
disease, cancer, blood, brain, …) with matching sample-info archives. A
versioned dataset registry, family-scoped downloads, explicit evaluation
metrics, and study-grouped holdouts are therefore scientific prerequisites
before trusting hierarchical comparisons.

## Decision

1. Insert **Milestone 5b** between Milestone 5 and Milestone 6:
   phenotype/source registry → family downloads → sample-info Parquet →
   evaluation metrics + study-grouped splits → TensorBoard on the flat loop →
   first multi-pack benchmark.
2. Amend ADR 0002 decision (4): the Stage 0 order is now
   … → flat baseline → **phenotype registry / multi-pack eval (5b)** →
   hierarchical → study-grouped cross-fitting.
3. Keep the Python package name `methyl-burden-score` and CLI entry point `mbs`
   unchanged. The public model name in docs and run metadata is **deepMAT**.
4. EWAS Atlas remains validation / enrichment only (ADR 0002 unchanged on
   Atlas vs Data Hub roles).

## Consequences

- [`docs/TODO_PIPELINE.md`](../TODO_PIPELINE.md) tracks Milestone 5b.
- Build brief: [`docs/plans/milestone-5b-phenotype-registry-eval.md`](../plans/milestone-5b-phenotype-registry-eval.md).
- Milestone 6 must not start until 5b acceptance criteria are met.
- DuckDB `phenotype` / `sample_phenotype` tables remain schema-ready but are
  not required to be fully populated in 5b (file registry is source of truth).
