# ADR 0010: Separate score export topology from phenotype comparators

## Status

Accepted

## Context

Milestone 7G named `C-mvalue-enet` the tissue phenotype winner, while the
intended product is a DeepRVAT-like reduced methylation representation:
sample×gene MBS plus non-gene features for downstream association. These are
different scientific outputs. Choosing the best tissue predictor must not
silently replace the score-export topology, and a useful score encoder need not
win every auxiliary phenotype.

The current 7F direct artifact is a fold-fitted contribution per task. The
association product instead requires direct CpG identity and values (or a
lossless reference into the canonical matrix).

## Decision

1. Final OOF product export uses the deepMAT cascade topology:
   gene-aggregated MBS, qualified orphan RBS kept per region, and indexed direct
   CpGs. TBS remains absent.
2. Phenotype benchmark tables rank methylation-input comparators separately.
   The current tissue comparator is `C-mvalue-enet`; the fold-selected
   `C-mvalue-enetS` benchmark is Milestone 7H.
3. Auxiliary phenotype heads train and select the shared encoder, but are not
   the downstream association product.
4. Claims must identify whether a result concerns product representation,
   phenotype prediction, or downstream association.
5. Final OOF cannot start until P4/P5 and 7H lock pooling, loss weights, fusion,
   direct-CpG export, and orphan-region eligibility.

## Orphan and direct semantics

- Orphan RBS are separate `region_id` columns, never one global score and never
  pooled by type.
- Only well-defined, versioned multi-CpG regions qualify.
- Region→gene allocation requires explicit evidence; unrestricted nearest-gene
  allocation is not accepted for final OOF.
- Singleton/unstructured loci remain direct.
- `direct_contrib.zarr` is a benchmark diagnostic; downstream association uses
  an indexed sample×direct-CpG block or lossless canonical-matrix view.

## Consequences

- A classical method may remain the best tissue comparator while deepMAT still
  produces the association representation.
- Reports cannot call a phenotype winner the “Milestone 7 topology.”
- Milestone 7H is a new pre-OOF gate.
- The final model remains array-manifest agnostic only in its shared
  aggregation path; task-specific direct weights are not agnostic to unseen
  loci.
