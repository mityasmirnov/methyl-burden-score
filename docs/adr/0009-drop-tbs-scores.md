# ADR 0009: Drop tile (TBS) scores; leftover CpGs stay direct

## Status

Accepted

## Context

[ADR 0006](0006-multipath-noncoding-scores.md) defined four product score
families: MBS, RBS, TBS, and direct CpG. Milestone 7C shipped graph-v2 with
50-CpG tile nodes; Milestone 7E trained independent TBS arms and fused
gene + RBS + TBS + direct as linear models on presence-aware region means.

Tiles randomly bin leftover loci into features that are not typed regulatory
regions. 7E evidence did not justify keeping TBS as a shipped score family.
Milestone **7F** locks leftover CpGs on the **direct** path instead.

Nearest-gene allocation of an already-typed **RBS** onto a gene (to form MBS)
is allowed. [ADR 0004](0004-unmapped-probe-retention.md) still forbids
collapsing an unmapped **CpG** into a nearest-gene proxy.

## Decision

1. **Product score families after 7F:**
   - **RBS** — one score per typed region (gene roles, CGI/shore, and later
     cCRE / enhancer / DMR / ChromHMM / similar);
   - **MBS** — gene-aggregated RBS (typed gene role and/or nearest-gene
     allocation of typed RBS);
   - **orphan RBS** — typed region scores with no gene allocation;
   - **direct** — per-locus contribution for CpGs with no typed-region
     assignment (Level-1 MAD z + fold-fitted elastic-net / sparse term).
2. **TBS scores are not a product arm.** Graph-v2 tile nodes may remain on
   disk unused. Train-time assignment ignores `region_system=tbs`.
3. Late fusion concatenates **saved** `[orphan RBS | MBS | direct]` matrices,
   not presence-aware region-mean tables.
4. ADR 0006’s assignment order is amended for scoring: annotation-first typed
   regions → RBS / MBS; leftover → **direct** (no adaptive tiles for scores).
5. ADR 0004 and the ban on nearest-gene for leftover CpGs stand.

## Consequences

- Milestone **7F** implements the cascade; **7G** and Milestone **7** export
  no TBS score matrices.
- [`DATA_CONTRACT.md`](../DATA_CONTRACT.md) marks `tbs.zarr` unused after 7F.
- [`ANNOTATION_GRAPH.md`](../ANNOTATION_GRAPH.md) documents score topology vs
  on-disk tile nodes.
- Frozen 7E runs and graph-v2 artifacts are not rewritten.

## Non-goals

- Deleting tile nodes from graph-v2.
- Retraining v0.1 or 7E checkpoints.
- Treating metadata-only predictors as methylation methods.
