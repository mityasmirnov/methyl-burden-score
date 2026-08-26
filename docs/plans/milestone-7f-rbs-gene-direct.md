# Plan: Milestone 7F — RBS→gene cascade + direct leftover

Status: **done** (implementation brief for Stage 0 Milestone **7F**).
Report: [`reports/inspection/stage0_7f_rbs_gene_direct/`](../../reports/inspection/stage0_7f_rbs_gene_direct/).
Normative: [ADR 0004](../adr/0004-unmapped-probe-retention.md),
[ADR 0006](../adr/0006-multipath-noncoding-scores.md),
[ADR 0009](../adr/0009-drop-tbs-scores.md),
[TODO_PIPELINE.md](../TODO_PIPELINE.md) §7F,
[post-v0-scientific-programme.md](post-v0-scientific-programme.md) §7F.

## Scope and acceptance

**Done when** (from TODO):

- Assignment + trainer implement the cascade on frozen `hub-ats-7e-3fold-v1`
- Fixture tests: leftover→direct; RBS→gene (nearest-gene allocation of RBS OK)
- Per-sample RBS / gene-RBS (MBS) / direct score matrices written and fused
  (not region-mean linear stand-ins)
- No TBS arm in the 7F model matrix
- Report: `reports/inspection/stage0_7f_rbs_gene_direct/`

**Budget:** prove topology + saved-score fusion. Hub smoke reuses the 7E
ceiling (2 epochs / 8192 loci). Raising that budget is **7G**.

## Locked decisions

| Choice | Decision | Why |
|--------|----------|-----|
| Graph | Keep `graph-grch38-gencode38-cgi-tile-v2`; no rebuild | Tile nodes unused; filter at train time |
| Typed regions | `region_system ∈ {gene, rbs}` | Gene five-role + CGI RBS; ignore TBS |
| RBS | One score per typed region | Glossary: gene roles are RBS too |
| MBS | Presence-aware max-pool of allocated RBS | Gene-aggregated RBS, not a separate gene FlatDeepSet |
| Nearest-gene | Same-chrom nearest gene for null-`gene_id` RBS; no gene on chrom → orphan | Allowed for typed RBS→MBS only |
| Direct | Leftover loci → Level-1 z + fold elastic-net | ADR 0006 transparent direct; no tiles |
| Fusion | Saved `[orphan_rbs \| mbs \| direct]` → linear heads | Not `presence_aware_means` |
| ADR | [0009](../adr/0009-drop-tbs-scores.md) drops product TBS | Amend 0006 product list |
| Splits | `hub-ats-7e-3fold-v1` only | Do not overwrite freezes |

## Schemas / contracts

Score export under `$MBS_ARTIFACT_ROOT/scores/<run_id>/`:

```text
mbs.zarr                  [n_samples, n_genes]
gene_present.zarr         [n_samples, n_genes]
rbs.zarr                  [n_samples, n_orphan_rbs]   # orphan only
direct_contrib.zarr       [n_samples, n_direct_tasks] # enet preds (age/tissue/sex)
sample_index.parquet
gene_index.parquet
region_index.parquet      # orphan RBS ids
score_manifest.json       # ADR 0008 polarity; no tbs.zarr
```

`direct_contrib` columns are fold-fitted elastic-net phenotype predictions on
leftover Level-1 z (one column per fitted task present in the fold). Fusion
concatenates orphan RBS + MBS + these columns.

Config: `configs/experiment/stage0_7f_rbs_gene_direct.yaml`.
CLI: `mbs train cascade`.

## Data / artifact flow

```text
graph-v2 + matrix-hub-age-tissue-sex-full-v1 + hub-ats-7e-3fold-v1
  → cascade_assign (exclude tbs; nearest-gene for null-gene RBS)
  → CascadeDeepSet train (CpG→RBS→MBS) + Level-1 + direct enet
  → write score Zarrs (no tbs)
  → late_fusion on saved matrices
  → reports/inspection/stage0_7f_rbs_gene_direct/
```

## Non-goals

- Full-budget bake-off / classical M-value comparators (**7G**)
- Final 5×6 OOF (**Milestone 7**)
- Graph rebuild / deleting tile nodes
- Shipping 7E region-mean multipath as the product topology

## Open questions

None (ADR 0009 closed the TBS product question).
