# Methylation Burden Score

A research codebase for learning a shared, gene-level DNA methylation burden score from variable sets of CpG measurements.

The public model name is **deepMAT** (deep Methylation Aggregation Transformer /
Deep Set family). The Python package remains `methyl-burden-score` with the
`mbs` CLI entry point; do not treat a package rename as required for Stage 0.

Stage 0 open training and pilot matrices use the CNCB **EWAS Data Hub** (with EWAS Atlas for later association checks). The model path is flat DeepRVAT-style **deepMAT** baseline → phenotype registry / multi-pack eval → real Hub pack matrices → multitask shared encoder (**5c done**) → **max-N flat age/tissue/sex (5d — done)** → hierarchical residual-path (**6 — done**) → **harmonized release + census (7A — current gate)** → nine-pack matrices / architecture / normalization / development CV (**7B–7E**) → study-grouped OOF cross-fitting (**7**, blocked until 7A–7E). Authoritative progress: [`docs/TODO_PIPELINE.md`](docs/TODO_PIPELINE.md). See also [`docs/STRATEGIC_PLAN.md`](docs/STRATEGIC_PLAN.md), [`docs/plans/post-v0-scientific-programme.md`](docs/plans/post-v0-scientific-programme.md), [`docs/adr/0002-ewas-datahub-primary-source.md`](docs/adr/0002-ewas-datahub-primary-source.md), [`docs/adr/0005-catalog-matrix-independence.md`](docs/adr/0005-catalog-matrix-independence.md), [`docs/adr/0006-multipath-noncoding-scores.md`](docs/adr/0006-multipath-noncoding-scores.md), and [`docs/adr/0007-crossfit-prerequisites.md`](docs/adr/0007-crossfit-prerequisites.md).

The Stage 0 implementation follows four design principles:

1. A DeepRVAT-style scoring function is shared across CpGs, regulatory regions, genes, and training traits.
2. The model consumes ragged CpG sets rather than a fixed array manifest.
3. CpGs are organized into biologically typed regions before gene-level aggregation.
4. Every reported training-sample score is obtained by study-grouped cross-fitting.

## Stage 0 scope

Stage 0 implements:

- canonical GRCh38 locus, probe, region, and gene registries;
- DuckDB/Parquet metadata catalogs (populated release = Milestone **7A**);
- Zarr matrix-store interfaces (backend-independent protocol in 7B);
- a flat CpG-to-gene Deep Set baseline;
- a hierarchical CpG-to-region-to-gene Deep Set model (v0.1 residual baseline);
- planned multi-path RBS / TBS / direct CpG scores (Milestone **7C**);
- static CpGPT sequence-adapter features exported offline;
- optional MethylGPT CpG-token priors as an ablation;
- age-regression and tissue-classification heads;
- study-grouped OOF cross-fitting after 7A–7E (Milestone **7**);
- array missingness and manifest-downsampling tests.

Stage 0 deliberately excludes dynamic foundation-model token extraction, LoRA, imputation during training, episignature classifiers, epivariant calling, production long-read training, ClickHouse, and default TileDB migration.

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

End-to-end CpG processing, flat vs hierarchical aggregation, phenotype-masked
training, and what “producing MBS” means today:
[`docs/SCORING_PIPELINE.md`](docs/SCORING_PIPELINE.md).

Probe / array annotation coverage (HM450, EPIC, EPICv2; assigned vs unassigned
shares): [`docs/PROBE_ANNOTATION_COVERAGE.md`](docs/PROBE_ANNOTATION_COVERAGE.md).

Datasets, on-disk sizes, sample counts, and trait harmonization:
[`docs/DATA_CATALOG.md`](docs/DATA_CATALOG.md).

Normative contracts: [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md),
[`docs/ANNOTATION_GRAPH.md`](docs/ANNOTATION_GRAPH.md).

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

Useful commands:

```bash
uv run mbs doctor --create-directories
uv run mbs catalog init
uv run mbs inspect ewas-metadata
uv run mbs train flat --overfit-fixture
# after Hub packs / GSE35069 are on disk under $MBS_DATA_ROOT:
# make download-ewas-study STUDY=GSE35069
# uv run mbs matrix convert --study-id GSE35069 --platform-id HM450 --verify
# uv run mbs matrix convert-pack --help
# uv run mbs train flat --config configs/experiment/stage0_flat_pilot.yaml
# Milestone 5c multitask (age + tissue):
# uv run mbs phenotypes build-multitask-table
# CUDA_VISIBLE_DEVICES=0 uv run mbs train flat \
#   --config configs/experiment/stage0_flat_multitask.yaml \
#   --run-id stage0-flat-multitask-age-tissue-v1
```

### Live monitoring (TensorBoard + TUI)

When `logging.tensorboard: true` (Hub / 5c configs), **`mbs train flat` starts
TensorBoard by default** (`logging.auto_tensorboard`, default on with TB). The
train JSON summary prints `tensorboard_url` and a `monitor_hint`.

In a **second SSH session** (TUI needs its own terminal):

```bash
source scripts/activate_data_environment.sh
uv run mbs monitor --run-id stage0-flat-multitask-age-tissue-v1
# starts/reuses TensorBoard + live Rich dashboard (epoch, loss, MAE, acc, GPU, ETA)
# --no-tensorboard   # TUI only
# --tb-port 6007     # if 6006 is taken by another run
```

Browser over SSH (local laptop):

```bash
ssh -L 6006:localhost:6006 <user>@<power-horse-host>
# open http://localhost:6006
```

If port 6006 is already in use, train/monitor pick the next free port and write
`$MBS_ARTIFACT_ROOT/runs/<run_id>/tensorboard.json` (URL + tunnel hint). Do not
start a second manual `tensorboard` for the same run.

Details: [`docs/EXPERIMENT_PROTOCOL.md`](docs/EXPERIMENT_PROTOCOL.md).

Hub/Atlas metadata contracts: [`docs/EWAS_METADATA.md`](docs/EWAS_METADATA.md).
Hub downloads: [`docs/EWAS_DATA.md`](docs/EWAS_DATA.md).
Data catalog (GB / samples / traits): [`docs/DATA_CATALOG.md`](docs/DATA_CATALOG.md).
Probe annotation coverage: [`docs/PROBE_ANNOTATION_COVERAGE.md`](docs/PROBE_ANNOTATION_COVERAGE.md).
Scoring pipeline schema: [`docs/SCORING_PIPELINE.md`](docs/SCORING_PIPELINE.md).
Inspection guide: [`docs/DATA_INSPECTION.md`](docs/DATA_INSPECTION.md).
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

Stage 0 milestones **1–6** are done (see
[`docs/TODO_PIPELINE.md`](docs/TODO_PIPELINE.md)):

| Done | What shipped |
|------|----------------|
| Annotation + static features | GRCh38 five-role graph; offline CpGPT locus features |
| Pilot matrix | GSE35069 EWAS_db → canonical Zarr (`mbs matrix convert`) |
| Flat deepMAT baseline | Overfit fixture + GSE35069 cell-type pilot train path |
| Phenotype registry (5b) | Versioned Hub packs, sample-info Parquet, study-grouped eval helpers |
| Hub metadata (5b′) | Atlas/Hub column contracts ([`docs/EWAS_METADATA.md`](docs/EWAS_METADATA.md)) |
| Real Hub packs (5b″) | `mbs matrix convert-pack`; study-holdout matrices + `stage0_hub_real_benchmark/` |
| Multitask deepMAT (5c) | Masked age+tissue heads on shared flat encoder |
| Max-N DeepRVAT flat (5d) | Uncapped age/tissue/sex GSM-union (`matrix-hub-age-tissue-sex-full-v1`, 13548 samples); run `stage0-flat-deeprvat-age-tissue-sex-full-v1`; report `reports/inspection/stage0_5d_max_n/` |
| Hierarchical residual path (6) | `mbs train hierarchical`; typed CpG→region→gene + residual slot (no `__unassigned__`); run `stage0-hier-deeprvat-age-tissue-sex-full-v1` (**deepMAT-hierarchical-v0.1**); report `reports/inspection/stage0_6_hierarchical/` |

**Current gate — Milestone 7A:** harmonized data release + phenotype census
([`docs/plans/post-v0-scientific-programme.md`](docs/plans/post-v0-scientific-programme.md)).
Final OOF cross-fitting (Milestone **7**) is blocked until 7A–7E. Milestone 6
closed: hierarchical vs flat on the same 5d folds; mapped path carries signal;
one-scalar residual_only is near chance for tissue/sex (bottleneck, not
biology). Flat 5d remains the stronger phenotype reference. Disease profile zip
is complete; EWAS_db All-Data mirror is still in progress (not required for
7A). See [`docs/DATA_CATALOG.md`](docs/DATA_CATALOG.md).

Public model name remains **deepMAT**; package/CLI stay `mbs` /
`methyl-burden-score`.

## Licensing

A project source-code license has not yet been selected. Reference papers, source repositories, datasets, and pretrained model weights have separate licenses and must be reviewed independently before redistribution or production use.