# ADR 0006: Multi-path gene, regulatory, tile, and direct CpG scores

## Status

Accepted

## Context

[ADR 0004](0004-unmapped-probe-retention.md) correctly requires retaining
unmapped CpGs without nearest-gene assignment. Milestone 6 implemented that as
a hierarchical **residual path** that max-pools all residual CpGs into **one
scalar** (~108k HM450 residual columns on the age/tissue/sex union). Residual-
only eval-time masking was near chance for tissue/sex; that tests a one-
dimensional bottleneck, not whether noncoding CpGs carry phenotype signal.

Flat deepMAT drops loci without a locus→region→gene path. Annotation coverage
shows ~22–30% of mapped probes/loci outside the five gene-region roles.
MethylCapsNet supports multiple capsule systems (genes, CGI, enhancers, bins);
its useful lesson is flexible region systems—not one catch-all leftover capsule.

DeepRVAT’s phenotype module is linear over **gene** impairment scores; a direct
CpG term is a deepMAT extension, not an exact DeepRVAT copy.

## Decision

1. **ADR 0004 stands** (retain unmapped; no nearest-gene). The Milestone 6 one-
   scalar residual path is frozen as **deepMAT-hierarchical-v0.1** only—it is
   not the preferred phenotype architecture going forward.
2. **Target score families** (Milestone 7C graph v2 + model):
   - **MBS** — gene methylation burden (`CpG → typed gene region → gene`);
   - **RBS** — non-gene regulatory burden (`CpG → cCRE / enhancer / CGI / DMR /
     ChromHMM / similar → regulatory region score`);
   - **TBS** — intergenic tile burden (adaptive CpG-count tiles for loci not
     assigned to gene or regulatory regions);
   - **Direct CpG** — per-locus (or context-generated) contribution for
     remaining / explicitly retained loci; no pre-compression of ~10⁵ loci into
     one scalar.
3. **Assignment order:** annotation-first (gene roles, then non-gene regulatory),
   then adaptive tiles, then direct CpG. Do not force every locus into a gene.
4. **Ablations** must train branches **independently** on identical folds.
   Eval-time masking of a jointly trained branch is insufficient evidence that
   a branch is uninformative.
5. Phenotype prediction may combine MBS + RBS + TBS + optional direct term +
   covariates. Prefer a transparent sparse per-locus direct branch as baseline;
   context-generated weights for cross-platform generalization.

## Consequences

- Graph release beyond five-role gene topology is Milestone **7C** (not §8
  “richer annotations” leftover). Product-surface cCRE portals remain deferred.
- [`ARCHITECTURE.md`](../ARCHITECTURE.md), [`ANNOTATION_GRAPH.md`](../ANNOTATION_GRAPH.md),
  and [`SCORING_PIPELINE.md`](../SCORING_PIPELINE.md) document today (v0.1) vs
  target (7C).
- Flat gene-only and hierarchical gene-only remain mandatory reference arms in
  Milestone **7E**.

## Non-goals

- Nearest-gene assignment for intergenic CpGs.
- One monolithic residual capsule as the production noncoding path.
- Requiring enhancer→gene evidence edges without a versioned evidence policy.
- Replacing or deleting the frozen hierarchical-v0.1 artifacts.
