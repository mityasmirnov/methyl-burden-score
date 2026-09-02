# ADR 0010: Typed-region → gene allocation policy

## Status

Accepted

## Context

Milestone 7F allows nearest-gene allocation of typed RBS with null `gene_id` onto
MBS. That is appropriate for product scoring sensitivity, but **7G′ Stage A**
architecture selection requires a fair gene-only benchmark: neural encoders and
`C-mvalue-*-G` must share CpGs that have **evidence-backed** region→gene edges,
not regions forcibly assigned to the nearest same-chromosome gene without a
distance or evidence threshold.

ADR 0004 forbids nearest-gene assignment for **unmapped CpGs** (leftover direct
path). This ADR governs **typed RBS → gene** allocation only.

## Decision

`build_cascade_assignment` accepts an explicit `gene_allocation` policy:

| Policy | Behaviour |
|--------|-----------|
| `explicit_only` | Keep annotation-backed `gene_id` only; null-gene typed regions stay orphan |
| `bounded_nearest` | Nearest same-chromosome gene only if distance ≤ `max_nearest_gene_bp` |
| `legacy_nearest` | Unconditional nearest-gene (7F default; sensitivity only) |

**7G′ Stage A** uses `explicit_only`. Stage B / product paths may use
`legacy_nearest` or `bounded_nearest` when documented.

Persist allocation policy and matched `gene_col_indices` in
`gene_panel_manifest.json` so neural and classical arms share identical columns.

## Consequences

- Stage A rerun uses new run IDs (e.g. `*-explicit`) — do not overwrite
  contaminated pre-fix checkpoints.
- Reports must refuse architecture lock until test-only `mbs_e2e` and completed
  classical `-G` folds exist on the same panel.
- Sensitivity analyses comparing `explicit_only` vs `legacy_nearest` quantify
  how much nearest-gene inflates the gene-only panel.

## Non-goals

- Do not revert ADR 0004 direct-path retention.
- Do not assign intergenic CpGs to nearest genes.
