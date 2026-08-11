# ADR 0004: Retain unmapped CpGs and probes in the hierarchical model

## Status
Accepted

## Context

Milestone 6 adds the hierarchical CpG→region→gene model. The repository already uses typed regulatory regions for annotated loci, but the next model stage must also preserve loci that are not mapped to a regulatory region.

The retention policy is:

- mapped loci keep their regulatory typing and gene assignment;
- unmapped CpGs/probes remain in the model;
- unmapped loci are not collapsed into the nearest gene;
- performance must be reported separately for mapped and unmapped loci;
- absence of annotation is not equivalent to absence of signal.

## Decision

The hierarchical model will keep all CpG/probe observations in the batch contract and model pipeline.

Annotated loci follow the typed hierarchy:

```text
CpG -> region role -> gene
```

Unmapped loci follow a separate retained path and are never silently merged into typed regions or nearest-gene proxies.

## Consequences

- Matrix conversion must preserve unmapped loci instead of filtering them away.
- The batch contract must carry annotation-status masks.
- The hierarchical model must expose a residual/unmapped path or equivalent retained output.
- Evaluation must stratify metrics by mapped versus unmapped loci.
- Milestone 6 planning and implementation must not use nearest-gene assignment for unannotated loci.

## Non-goals

- Do not introduce a new nearest-gene intergenic policy.
- Do not require dynamic foundation-model token extraction for this stage.
- Do not treat unannotated loci as absent or noise by default.
