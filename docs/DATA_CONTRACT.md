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

```text
sample_id
study_id
source_sample_id
donor_id
replicate_group
age
sex
tissue_raw
tissue_ontology_id
case_control
metadata_json
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

Large sample×locus arrays live in **Zarr/HDF5**, not DuckDB. DuckDB holds
catalog metadata (and optional views over Parquet); phenotype labels for
training are Parquet (`canonical/phenotypes/`). See Milestone 5c storage notes
in [`plans/milestone-5c-multitask-shared-encoder.md`](plans/milestone-5c-multitask-shared-encoder.md).

Each canonical study matrix has:

```text
betas.zarr or betas.h5       [n_samples, n_study_loci]
sample_index.parquet         study row -> sample_id
locus_index.parquet          study column -> canonical locus_id
matrix_manifest.json         provenance and array metadata
```

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
mbs.zarr                  [n_samples, n_genes]
gene_present.zarr         [n_samples, n_genes]
gene_observed_count.zarr  [n_samples, n_genes]
region_scores.zarr        optional
sample_index.parquet
gene_index.parquet
score_manifest.json
```

The score manifest records fold and restart membership. Every sample must be traceable to models that excluded its study and replicate group.

## Inspection policy

Agents and routine code review use aggregate reports under `reports/inspection/`.
Do not embed raw methylation matrices in docs. Hub sample-info `.txt` extracts
under `reports/inspection/ewas_datahub_samples/` and Atlas small TSVs are
allowed for phenotype contracts; follow [`DATA_INSPECTION.md`](DATA_INSPECTION.md)
and [`EWAS_METADATA.md`](EWAS_METADATA.md) (`mbs inspect ewas-metadata`).
