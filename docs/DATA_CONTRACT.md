# Data contract

## Design goals

The canonical representation must support:

- 450K, EPIC, and EPICv2 arrays;
- variable observed CpG sets;
- study-grouped cross-validation;
- static locus-feature lookup;
- many-to-many locus-to-region and region-to-gene edges;
- future long-read coverage and count fields;
- reproducible conversion from immutable raw sources.

## Stable identifiers

### Sample

```text
sample_id      project-controlled stable identifier
source_id      source accession or supplied identifier
donor_id       biological donor when known
replicate_id   technical or longitudinal replicate group when known
study_id       canonical study or cohort identifier
```

Do not use a row number as the only stable sample identifier.

### Locus

Canonical locus key:

```text
GRCh38:chr{chromosome}:{cytosine_position}
```

The manifest records:

- genome build;
- one-based or zero-based convention;
- cytosine rather than probe-start position;
- strand normalization policy;
- source probe IDs;
- liftover provenance;
- ambiguity status.

### Region and gene

```text
region_id = stable graph-release identifier
gene_id   = stable Ensembl gene identifier without version suffix
gene_name = display symbol only; never the primary key
```

## Catalog tables

DuckDB holds metadata and may query Parquet. Large sample×locus arrays stay in
Zarr ([ADR 0005](adr/0005-catalog-matrix-independence.md)). Milestone **7A**
populates these tables into a versioned release via
`mbs catalog refresh-release`. `mbs catalog init` remains schema-only at
`canonical/catalog/`. See
[`plans/milestone-7a-harmonized-release.md`](plans/milestone-7a-harmonized-release.md).

### Versioned release layout (7A)

```text
$MBS_DATA_ROOT/canonical/releases/deepmat-data-v1/
├── release_manifest.json
├── catalog/
│   ├── catalog.duckdb
│   └── tables/
├── matrices/
├── phenotypes/
├── ontologies/
├── annotations/
├── graphs/
└── splits/
```

`release_manifest.json` records source checksums, retrieval dates, preprocessing
level, probe universe, genome build, graph/static-feature versions, phenotype
families, sample-dedup decisions, and code commit.

### `study`

```text
study_id
source_release_id
gse_id
cohort_id
platform_id
processing_level
genome_build
retrieved_at
```

### `sample`

Identity and study linkage. Phenotype values are **not** the long-form SoT
(optional denormalized age/sex/tissue_raw may exist for convenience):

```text
sample_id
study_id
source_sample_id
donor_id
replicate_group
age                 # optional denormalized convenience
sex                 # optional
tissue_raw          # optional
tissue_ontology_id  # optional
case_control        # optional
metadata_json
```

### `sample_source_membership` (7A)

One GSM may appear in several Hub packs. One biological `sample_id`, many
memberships:

```text
sample_id
source_release_id
phenotype_family
source_file
source_row
matrix_id
row_index
```

### `assay_file`

```text
assay_file_id
study_id
path
format
sha256
n_rows
n_columns
matrix_orientation
schema_hash
processing_level
```

### `phenotype`

```text
phenotype_id
phenotype_name
phenotype_type
ontology_id
unit
```

### `sample_phenotype` (long-form SoT)

Do not expand one wide row with hundreds of disease columns. Multi-label and
multi-pack observations use long form; include `source_family` so duplicates
are not overwritten. Missing disease labels are **unknown**, not automatic
controls, unless pack documentation establishes controls.

```text
sample_id
phenotype_id
numeric_value          # SQL column name (not value_numeric)
categorical_value      # SQL column name (not value_categorical)
label_status           # observed | unknown | control | case
is_observed
source_family
source_record_id
ontology_id
```

Observation uniqueness: `(sample_id, phenotype_id, source_family)` or an
explicit `observation_id`. Wide `sample_phenotype_table.parquet` used for
training remains a **derived** join table.

### `probe`

```text
probe_id
platform_id
probe_design
manifest_version
quality_flags
```

### `locus`

```text
locus_id
genome_build
chromosome
position
canonical_key
mapping_status
```

### `probe_locus_edge`

```text
probe_id
locus_id
mapping_source
mapping_confidence
is_primary
```

### `region`

```text
region_id
region_type
chromosome
start
end
strand
gene_id
source
source_version
```

### `locus_region_edge`

```text
locus_id
region_id
edge_weight
evidence_type
primary_gene_role
```

## Numeric matrix contract

Large sample×locus arrays live in **Zarr** (or a future `MethylationStore`
backend), not DuckDB. DuckDB holds catalog metadata and optional views over
Parquet; phenotype labels for training are Parquet (`canonical/phenotypes/`).
See [ADR 0005](adr/0005-catalog-matrix-independence.md) and Milestone 5c notes
in [`plans/milestone-5c-multitask-shared-encoder.md`](plans/milestone-5c-multitask-shared-encoder.md).

### Matrix store (Zarr layout; protocol deferred)

Training currently opens the Zarr layout via `mbs.matrix.store`. A thin
`MethylationStore` protocol remains deferred (YAGNI until a second backend is
needed; [ADR 0005](adr/0005-catalog-matrix-independence.md)). TileDB / Parquet
observation stores are future backends only.

Each canonical study matrix has:

```text
betas.zarr                   [n_samples, n_study_loci]  (compressed float32)
sample_index.parquet         study row -> sample_id
locus_index.parquet          study column -> canonical locus_id
sample_phenotypes.parquet    long-form labels (may repeat sample_id)
matrix_manifest.json         provenance and array metadata
```

`locus_index.parquet` may include Milestone **7B** columns:

```text
contributing_probe_ids       pipe-separated probe IDs collapsed into the column
collapse_method              identity | mean | median
```

### Virtual Hub pack index (7B)

Cross-pack GSM membership without a dense nine-pack union:

```text
canonical/matrices/hub_pack_matrix_index.parquet
  family, matrix_id, sample_id, row_index, platform, betas_path
```

Overlapping GSM betas are **verified** across packs (max abs diff on a locus
subset). Do not silently take the first pack when merging; discordant pairs
block a claimed merged beta. See
[`plans/milestone-7b-complete-hub-matrices.md`](plans/milestone-7b-complete-hub-matrices.md).

The matrix manifest contains:

```text
artifact_version
source_file_hashes
conversion_commit
study_id
platform_id
processing_level
genome_build
dtype
shape
chunking
compression
missing_value_encoding
sample_index_hash
locus_index_hash
created_at
```

## Methylation values

Raw beta values are retained for Stage 0.

Validation rules:

- finite observed beta values normally lie in `[0,1]`;
- out-of-range values are reported, not silently clipped;
- missing values remain explicit;
- M-values are derived lazily or stored in a separate artifact;
- train-fold standardization parameters never modify canonical raw beta values.

M-value transformation:

```math
M=\log_2\frac{\beta+\epsilon}{1-\beta+\epsilon}
```

The epsilon value is part of the resolved experiment configuration.

## Static-feature contract

Static features are stored once per canonical locus:

```text
embeddings.zarr or embeddings.npy    [n_loci, d]
loci.parquet                          embedding row -> locus_id
artifact.json                         feature provenance
```

The manifest contains:

```text
feature_set_id
source_model
source_repository
source_commit
checkpoint_sha256
configuration_sha256
vocabulary_sha256, when applicable
source_locus_table_sha256
genome_build
input_dimension
output_dimension
storage_dtype
normalization
export_command
created_at
```

A training job fails if the feature-locus hash differs from the active locus registry.

## Graph contract

A graph release contains:

```text
genes.parquet
regions.parquet
locus_region_edges.parquet
region_gene_edges.parquet
graph_manifest.json
regions.bed
```

The manifest records all source annotation versions, interval policies, overlap precedence, genome build, input hashes, and graph-builder commit.

## Ragged batch contract

```python
@dataclass
class MethylationBatch:
    sample_ids: list[str]
    beta: Tensor
    locus_row: Tensor
    sample_index: Tensor
    edge_cpg_index: Tensor
    edge_region_index: Tensor
    region_type: Tensor
    region_gene_index: Tensor
    gene_sample_index: Tensor
    gene_panel_index: Tensor
    targets: dict[str, Tensor]
    target_masks: dict[str, Tensor]
```

The batch stores observed sample–CpG pairs, not a padded whole-manifest matrix.

## Score artifact contract

Out-of-fold scores are stored as:

```text
mbs.zarr                  [n_samples, n_genes]   # gene-aggregated RBS (7F+)
gene_present.zarr         [n_samples, n_genes]
gene_observed_count.zarr  [n_samples, n_genes]
rbs.zarr                  optional [n_samples, n_orphan_or_all_rbs]
tbs.zarr                  unused after 7F (ADR 0009); may be absent
direct_contrib.zarr       optional per-task direct CpG contributions
sample_index.parquet
gene_index.parquet
score_manifest.json
```

From Milestone **7F**, product fusion uses saved **orphan RBS + MBS + direct**
matrices ([ADR 0009](adr/0009-drop-tbs-scores.md)); do not write a TBS score
arm. The score manifest records fold and restart membership, **score polarity /
orientation anchor** ([ADR 0008](adr/0008-score-identifiability.md)), and
fold-fitted normalizer hashes. Every sample must be traceable to models that
excluded its study and replicate group. Milestone **7** (after 7F and 7G) is
the OOF export gate. Averaging unaligned `MBS` and `1−MBS` folds is undefined.
OOF MBS is a **predictive** sample×gene representation, not a constraint or
LOEUF analogue.

## Trait eligibility (7A)

Derived Parquet / DuckDB table `trait_eligibility` for every harmonized
phenotype (columns: `phenotype_id`, family, task type, n_samples / cases /
controls, prevalence, n_studies / platforms / tissues, confounding scores,
`eligible_core_task` / `eligible_auxiliary_task` / `eligible_external_evaluation`,
`exclusion_reason`). Cutoffs:
[`plans/post-v0-scientific-programme.md`](plans/post-v0-scientific-programme.md).

## Inspection policy

Agents and routine code review use aggregate reports under `reports/inspection/`.
Do not embed raw methylation matrices in docs. Hub sample-info `.txt` extracts
under `reports/inspection/ewas_datahub_samples/` and Atlas small TSVs are
allowed for phenotype contracts; follow [`DATA_INSPECTION.md`](DATA_INSPECTION.md)
and [`EWAS_METADATA.md`](EWAS_METADATA.md) (`mbs inspect ewas-metadata`).
