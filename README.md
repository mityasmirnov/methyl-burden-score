# Methylation Burden Score

A research codebase for learning a shared, gene-level DNA methylation burden score from variable sets of CpG measurements.

The Stage 0 implementation follows four design principles:

1. A DeepRVAT-style scoring function is shared across CpGs, regulatory regions, genes, and training traits.
2. The model consumes ragged CpG sets rather than a fixed array manifest.
3. CpGs are organized into biologically typed regions before gene-level aggregation.
4. Every reported training-sample score is obtained by study-grouped cross-fitting.

## Stage 0 scope

Stage 0 implements:

- canonical GRCh38 locus, probe, region, and gene registries;
- DuckDB/Parquet metadata and annotation catalogs;
- Zarr or HDF5 matrix-store interfaces;
- a flat CpG-to-gene Deep Set baseline;
- a hierarchical CpG-to-region-to-gene Deep Set model;
- static CpGPT sequence-adapter features exported offline;
- optional MethylGPT CpG-token priors as an ablation;
- age-regression and tissue-classification heads;
- study-grouped cross-fitting and leakage controls;
- array missingness and manifest-downsampling tests.

Stage 0 deliberately excludes dynamic foundation-model token extraction, LoRA, imputation during training, episignature classifiers, epivariant calling, and production long-read training.

## Server layout

All durable and transient files must remain under `/data`. Defaults are project-local so bootstrap does not need shared `/data/datasets` ownership:

```text
$MBS_ROOT                 Git working tree
$MBS_DATA_ROOT            data/ (canonical + staging)
$MBS_SCRATCH_ROOT         scratch/ (temporary computation)
$MBS_CACHE_ROOT           cache/ (project + tool caches)
$MBS_ARTIFACT_ROOT        artifacts/ (runs, checkpoints, scores)
```

Do not place datasets, environments, checkpoints, model weights, or caches under `$HOME`.

## Quick start on `power-horse`

```bash
mkdir -p /data/projects
cd /data/projects
git clone git@github.com:mityasmirnov/methyl-burden-score.git
cd methyl-burden-score

cp .env.example .env
source scripts/activate_data_environment.sh

uv sync --all-groups
uv run mbs doctor
uv run pytest
```

The repository does not contain research data, pretrained weights, or copies of the reference repositories. See [`docs/WORKSPACE.md`](docs/WORKSPACE.md) and [`scripts/add_reference_submodules.sh`](scripts/add_reference_submodules.sh).

## Architecture

For sample `s`, CpG `c`, regulatory region `r`, and gene `g`:

```text
CpG features
    -> shared CpG encoder
    -> permutation-invariant region pooling
    -> shared region encoder
    -> permutation-invariant gene pooling
    -> shared scalar compression network
    -> MBS[s, g]
```

The exact DeepRVAT-compatible flat baseline is also retained:

```text
CpG -> shared phi -> elementwise max by gene -> shared rho -> sigmoid MBS
```

Phenotype heads operate on the complete vector of gene scores and are not part of the exported MBS scoring function.

## Development commands

```bash
make doctor
make lint
make typecheck
make test-fast
make test
make catalog-init
make catalog-build
```

Useful first-pass commands after bootstrap:

```bash
uv run mbs doctor --create-directories
uv run mbs catalog init
uv run mbs inspect source --source-id cpgcorpus
uv run mbs inspect cpgcorpus-gpl --gse GSE125367 --gpl GPL21145
```
Machine-specific paths belong in `.env` or `configs/local/`; neither is committed.

## Repository policy

Committed:

- Python source;
- SQL schemas and views;
- YAML configurations;
- documentation and architecture decisions;
- small synthetic test fixtures;
- artifact manifests and checksums.

Never committed:

- sample-level methylation data;
- IDAT, BAM, VCF, Arrow, Parquet, Zarr, HDF5, or SQLite data artifacts;
- pretrained checkpoints or embeddings;
- secrets and credentials;
- generated run outputs;
- Docker layer data.

## Status

The repository is in Stage 0 bootstrap. Interfaces and tests are being established before full data ingestion and model training.

## Licensing

A project source-code license has not yet been selected. Reference papers, source repositories, datasets, and pretrained model weights have separate licenses and must be reviewed independently before redistribution or production use.