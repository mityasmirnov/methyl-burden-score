# ADR 0012: Seed-gene discovery vs deployment CpG input

## Status

Accepted

## Context

Fold-safe association for `internal_fold` seed panels uses a univariate
prefilter (today 4,096 CpGs) plus stability selection. Reports and informal
readouts often treat `n_seed_cpgs: 4096` as the deepMAT input width. That
misstates the contract: those CpGs **discover seed genes**; the G2/C2 arms and
the production encoder consume **all explicit gene-linked CpGs** of the selected
genes (siblings included), and deployment must accept whatever eligible CpGs a
platform observes for each gene.

ADR 0011 defines *where* seed genes come from (`external_clean` /
`internal_fold` / `hybrid_fold`). This ADR defines what the discovery CpG set
is for, how expanded gene CpGs relate to it, and what the deployment input
contract is. Parent plan:
[`docs/plans/milestone-7g-prime-age-seed-mask.md`](../plans/milestone-7g-prime-age-seed-mask.md).

## Decision

1. **Discovery ≠ input panel.** Fold-safe CpG association (prefilter +
   stability) ranks genes. The discovery-set size (e.g. 4,096) must never be
   reported as the number of CpGs supplied to G2, C2, or a production encoder.
2. **Gene enrichment is required for matched seed arms.** After seed genes are
   chosen, include every ADR 0010 `explicit_only` gene-linked CpG for those
   genes. Persist `is_seed_cpg` so discovery vs sibling loci remain auditable.
3. **Required panel fields (per trait)** — reports must distinguish:

   | Field | Meaning |
   |-------|---------|
   | `n_discovery_cpgs` | CpGs surviving stability selection |
   | `n_seed_genes` | selected genes |
   | `n_expanded_gene_cpg_edges` | all selected gene–CpG edges |
   | `n_unique_expanded_gene_cpgs` | unique CpGs used by G2/C2 |
   | `n_multigene_cpgs` | CpGs attached to multiple selected genes |
   | `seed_fraction_of_expanded` | unique discovery CpGs / unique expanded |

   Legacy aliases (`n_seed_cpgs`, `n_seed_cpgs_after_stability`) may remain for
   back-compat but must not be the primary human readout.
4. **Deployment contract.** Production deepMAT aggregates whatever eligible
   CpGs are observed for a gene (genomic coordinate + build/version, not probe
   ID). Pool only observed loci; never treat absent CpGs as methylation zero.
   Presence-aware / coverage-aware pooling is the reuse path
   (`presence_aware_means`, observed-edge feature assembly). Platform-specific
   transforms (array M-value vs ONT frequency) and coverage dropout robustness
   experiments are required before claiming platform-agnostic behaviour — they
   are **not** implemented by this ADR.
5. **Trait list is config-driven.** Seed-panel traits come from experiment
   YAML, not a hard-coded triple. ATS remains age / tissue / sex (+ optional
   `sex_autosome` control). New traits require labels on the training matrix,
   eligibility under the post-v0 census bar, and a head — not a rewrite of gene
   enrichment.

## Consequences

- G2 = expanded gene CpGs + trait seed masks; C2 uses the same expanded CpGs.
- Audit scripts and analysis tables must show discovery vs expanded side by
  side so `4096` cannot be read as G2 width.
- Adding BMI / disease / cancer is a later config + labels + head change;
  unknown disease/cancer labels stay unknown (not controls).

## Non-goals

- Implementing 450K↔EPIC dropout augmentation, ONT scale transforms, or
  coverage-confidence export tensors.
- Genericizing `MultitaskHeads` beyond age / tissue / sex.
- Changing ADR 0011 seed sources or ADR 0010 gene allocation.
