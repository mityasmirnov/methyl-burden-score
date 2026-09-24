# ADR 0011: Seed gene sources for phenotype-masked heads

## Status

Accepted

## Context

Milestone 7G introduces optional **seed masks** on the phenotype heads
(`MultitaskHeads`: age, tissue, sex) so a trait head reads only a curated subset
of gene columns rather than the full gene panel. A seed mask is a leakage-bearing
choice: which genes a head is *allowed* to use is itself learned or imported
information, and the scientific invariants in `AGENTS.md` (§2, §3) require any
phenotype-derived feature selection to be fitted **inside the training fold**.

We therefore need an explicit, auditable taxonomy of where seed genes come from
and which of those sources are fold-safe. This ADR fixes that taxonomy. The
parent build plan is
[`docs/plans/milestone-7g-prime-age-seed-mask.md`](../plans/milestone-7g-prime-age-seed-mask.md).

## Decision

Seed masks / gene panels are constructed from one of the sources below,
recorded per arm and per fold (or per catalog artifact):

| Source | Construction | Leakage |
|--------|--------------|---------|
| `external_clean` | Atlas panel after excluding overlapping benchmark studies/samples | Fixed external prior (**not** fold-fitted) |
| `internal_fold` | Associations / elastic-net importance computed on the outer **training fold only** | Fold-fitted |
| `hybrid_fold` | Combine `external_clean` + `internal_fold` **inside the training fold** | Fold-safe when combined inside train |
| `catalog_universe` | All-data univariate EWAS (± Atlas union) to build a **working probe/gene universe** for inspection and as a pre-filter before fold reselect | **Not** a CV training panel — uses labeled samples across folds |

Rules that bind these sources:

1. **Atlas genes are NOT automatically fold-fitted.** `external_clean` is a fixed
   external prior; treating it as fold-fitted would misreport the leakage class.
   It is only clean because overlapping benchmark studies/samples were excluded
   when the Atlas panel was derived — that exclusion is the whole point of the
   `_clean` suffix.
2. **Do not trust Atlas gene symbols as canonical allocation.** Atlas ships gene
   *symbols*; the pipeline allocates on annotation-graph gene IDs. Every seed
   gene must be remapped through the annotation graph using ADR 0010
   `explicit_only` allocation. Symbols with no `explicit_only` edge are dropped,
   not force-mapped to a nearest gene.
3. `internal_fold` and the `internal_fold` component of `hybrid_fold` must be
   derived using only outer-training-fold samples. No validation or test sample
   may inform the mask.
4. The resolved mask (source, gene count, and a content hash) is persisted in the
   run metrics / `score_manifest` so the leakage class of every reported number
   is auditable.
5. **`catalog_universe` is never the training feature set for CV.** Stage A/B
   (all-data catalog → restrict to genes associated with any trait) may shrink
   the probe universe so fold EWAS does not scan orphans / unassociated genes.
   Stage C must **reselect on outer-train only** (`internal_fold`) and expand
   seeds → gene-linked sibling CpGs (ADR 0012). Reports must label catalog
   artifacts `leakage_class: catalog_universe`.

## Consequences

- Reports must state the seed source per arm; an `external_clean` arm may not be
  described as fold-fitted, and an `internal_fold` / `hybrid_fold` arm may not be
  described as a fixed external prior. A `catalog_universe` artifact may not be
  described as the CV panel.
- A mask that selects too few genes is a configuration error, not a valid tiny
  panel: `MultitaskHeads` fails closed when a provided mask selects fewer than 32
  genes for any trait.
- Atlas-derived masks require an `explicit_only` remap step before use; the
  count of dropped (unmappable) symbols should be reported.
- Milestone **11** trait-universe runners
  (`scripts/run_trait_universe_catalog.py`,
  `scripts/run_trait_universe_fold_reselect.py`) implement catalog → restrict →
  fold-safe reselect; see
  [`plans/milestone-11-fold-selected-panel.md`](../plans/milestone-11-fold-selected-panel.md).

## Non-goals

- This ADR does not change the leftover direct-CpG path (ADR 0004) or the
  typed-region → gene allocation policy (ADR 0010); it consumes the latter.
- It does not mandate seed masks for every arm — dense (all-ones) heads remain
  the default when no mask is supplied.
