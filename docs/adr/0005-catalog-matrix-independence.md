# ADR 0005: Catalog independent of matrix store; DuckDB + Parquet + Zarr

## Status

Accepted

## Context

Stage 0 has a DuckDB schema, Parquet phenotype tables, Zarr matrices, and
inspection reports, but `mbs catalog init` / `build` only apply SQL. There is
no end-to-end job that populates catalog tables from the current Hub assets.
Sources of truth are fragmented across YAML, Parquet, Zarr manifests, and
reports. Nine Hub packs are downloaded; only age/tissue/sex (and their GSM
union) are converted to full canonical matrices.

Large sample×CpG arrays must not become long DuckDB tables. Future WGBS may
need a different array backend than Zarr, without rewriting the phenotype or
provenance catalog.

## Decision

1. **Keep three layers:** DuckDB (metadata, provenance, phenotypes, splits,
   artifact inventory), Parquet (durable tabular SoT), Zarr (dense numerical
   matrices). Do **not** store tens of billions of sample–CpG observations in
   DuckDB.
2. **Catalog ⊥ matrix store.** Training and conversion depend on a thin
   `MethylationStore` protocol (to be introduced when Milestone 7B touches the
   store). First implementation wraps the existing Zarr layout. TileDB sparse
   (or sharded Zarr v3) is benchmarked only when a representative WGBS cohort
   exists. **Do not add ClickHouse** for Stage 0.
3. **Versioned data releases** under
   `$MBS_DATA_ROOT/canonical/releases/deepmat-data-v1/` with
   `release_manifest.json` (source checksums, retrieval dates, preprocessing,
   probe universe, genome build, graph/static versions, phenotype families,
   dedup decisions, code commit) and a populated `catalog.duckdb`.
4. **Populate** catalog tables from Parquet (planned CLI:
   `mbs catalog refresh-release` / `validate-release`). Long-form
   `sample_phenotype` is the phenotype SoT; include `source_family` (and
   related observation keys) so multi-label / multi-pack rows are not
   overwritten. Add `sample_source_membership` for one GSM in several Hub
   packs. Wide training join tables remain *derived* views, not the SoT.
5. Identity fields stay on `sample`; denormalized `age` / `sex` /
   `tissue_raw` on `sample` are optional convenience only.

## Consequences

- Milestone **7A** delivers the first populated release and phenotype census
  ([`TODO_PIPELINE.md`](../TODO_PIPELINE.md);
  [`plans/post-v0-scientific-programme.md`](../plans/post-v0-scientific-programme.md)).
- [`DATA_CONTRACT.md`](../DATA_CONTRACT.md) sketches release layout and new
  tables; SQL/schema implementation lands with 7A.
- Frozen v0 assets (`deepmat-data-age-tissue-sex-v1` and flat/hier v0.1 runs)
  must not be overwritten by the new release job.

## Non-goals

- Migrating matrices to TileDB or ClickHouse now.
- Materializing sample×CpG long tables in DuckDB.
- Advertising `refresh-release` as implemented before the CLI exists.
