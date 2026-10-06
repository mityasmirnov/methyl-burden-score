# Milestone 12b §3.5 — sparse gather + epoch within-gene sampler

> Status (2026-10-06): **in progress**. K=16 fixed-cap smoke validated plumbing
> only — **not** §3.5 done. Product path = sparse/ragged memory + stochastic
> within-gene sampling. Parent:
> [`milestone-12b-full-gene-panel.md`](milestone-12b-full-gene-panel.md) ·
> seed recipe: [`milestone-12b-deeprvat-seed-recipe.md`](milestone-12b-deeprvat-seed-recipe.md).

## Scope and acceptance

**Done when:**

1. N-light (and cascade) training paths do **not** materialize
   `[n_samples × n_universe]` (`betas_ram = betas[:, :n_cols]` gone for
   flat_region product path).
2. Minibatch gather from `RoutedBetas`: sample rows × unique locus columns of
   the **current** sampled gene graph; nine-pack routing preserved; IO
   instrumentation (logical cols, physical cols/chunks, bytes, IO time,
   packing time).
3. Within-gene sampling runs **before** auto batch-size / VRAM calibration;
   calibrator probes the sampled graph. “440k→261k edges but batch still 55”
   is **not** an accepted end state.
4. Epoch-varying sampler: never drop a selected gene; ≤K use all; >K resample
   without replacement by epoch; same locus set reusable across samples in that
   epoch; persist exposure (`times_available`, `times_selected`,
   `selection_probability`). Val + score **never** sample.
5. Full-CpG scoring: train may use K; inference uses all observed linked CpGs;
   emit MBS, present, `n_observed_cpg`, `n_available_cpg`.
6. K∈{8,12,16,32} one-fold benchmark ranked on nested-enet + runtime (not
   single-seed `mbs_e2e`).
7. Product default after gather works: max samples · CpGPT on · multi-trait
   **seed-bank union** · trait heads; train φ/ρ on seed genes; score ~20k.
8. All-gene full-width train = **ablation**, not automatic production recipe.
9. Do **not** launch cascade 5×6 or dense full-width converge until above lands.

## Locked decisions

| Choice | Decision | Why |
|--------|----------|-----|
| Product train genes | Multi-trait seed union (DeepRVAT) | G1 PASS → encoder transfers; cheaper than 20k-gene train |
| Full-width all-gene train | Ablation only | Ask “does it beat seed-bank?” once gather works |
| Cap plumbing (`cap_flat_region_cpgs_per_gene`) | Keep as helper; train uses epoch sampler | Fixed subset starves large genes |
| Val / score | Full CpGs, never sample | Deterministic product export |
| `RoutedBetas.wide` | Contiguous-prefix only (no half-universe rule) | Sparse K graphs must not read full pack rows |
| Cascade 5×6 / dense converge | Blocked | Recipe incomplete |

## Schemas / config

```yaml
training:
  max_cpgs_per_gene: 12          # epoch-varying when set
  sparse_betas: true             # no betas_ram preload (flat_region default on)
  # gene_holdout.method: seed    # product after gather smoke
```

Exposure artifact (under run `scores/` or `reports/`): JSON/parquet with per
locus_id / edge: `times_available`, `times_selected`, `selection_probability`.

## Data / artifact flow

```text
TRAIN (product)
  seed genes (outer-fold) → full linked edges → epoch sample ≤K/gene
       → batch rows × unique sampled cols from RoutedBetas
       → φ/ρ + trait heads

SCORE
  all genes × ALL observed linked CpGs → frozen φ/ρ → MBS + present + coverage
```

## Non-goals / deferred

- Cascade 5×6, G3 platform dropout, ONT adapter, role-stratified sampling
  (later), gene-only-only product (keep MBS/RBS/direct paths).

## Open questions

- Optimal K after gather (resolve via §6 benchmark).
- Whether all-gene train ablation beats seed-bank (resolve after gather).
