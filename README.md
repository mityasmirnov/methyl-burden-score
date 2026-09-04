# Methylation Burden Score

A research codebase for learning a shared, gene-level DNA methylation burden
score (MBS) from variable sets of CpG measurements.

The public model name is **deepMAT** (deep Methylation Aggregation Transformer /
Deep Set family). The Python package remains `methyl-burden-score` with the
`mbs` CLI entry point.

Primary open data source: CNCB **EWAS Data Hub** (EWAS Atlas for association
checks). Authoritative progress:
[`docs/TODO_PIPELINE.md`](docs/TODO_PIPELINE.md).

**Current gate (2026-09-04):** **7G′ Stage B GPU** — fold-selected panel + full
model on the locked Stage A gene encoder. Stage A DeepRVAT Tier-1 screen is
**done**; provisional lock is **`P2-G`** (CascadeDeepSet, max/max pooling, 15
epochs) on the `explicit_only` gene-linked panel (51 375 CpGs). Cascade is
**not** ≥0.03 tissue F1 ahead of classical `C-mvalue-enet-G` (0.388 vs 0.373
`mbs_e2e`). Final Milestone **7** 5×6 OOF remains blocked until Stage B lands.

Programme docs: [`docs/STRATEGIC_PLAN.md`](docs/STRATEGIC_PLAN.md),
[`docs/plans/post-v0-scientific-programme.md`](docs/plans/post-v0-scientific-programme.md),
[`docs/plans/milestone-7g-prime-matched-probe-lightweight.md`](docs/plans/milestone-7g-prime-matched-probe-lightweight.md),
[`docs/plans/milestone-7g-prime-stage-a-deeprvat-screen.md`](docs/plans/milestone-7g-prime-stage-a-deeprvat-screen.md).
ADRs: [0002](docs/adr/0002-ewas-datahub-primary-source.md) (Hub primary),
[0007](docs/adr/0007-crossfit-prerequisites.md) (OOF),
[0008](docs/adr/0008-score-identifiability.md) (orientation),
[0009](docs/adr/0009-drop-tbs-scores.md) (no TBS),
[0010](docs/adr/0010-gene-allocation-policy.md) (`explicit_only`).

Do **not** retrain frozen **deepMAT-flat-v0.1** / **hierarchical-v0.1**.

## Design principles

1. A DeepRVAT-style scoring function is shared across CpGs, typed regions,
   genes, and training traits.
2. The model consumes ragged CpG sets rather than a fixed array manifest.
3. CpGs are organized into biologically typed regions before gene-level
   aggregation (**RBS → gene MBS**); leftover CpGs stay **direct** ([ADR 0009](docs/adr/0009-drop-tbs-scores.md) — **no tile/TBS scores**).
4. Every reported training-sample score is obtained by study-grouped
   cross-fitting (Milestone **7**, after 7G′).

## Stage 0 scope

Stage 0 implements:

- canonical GRCh38 locus, probe, region, and gene registries (graph-v2 on disk);
- DuckDB/Parquet metadata catalogs (populated release = Milestone **7A**);
- Zarr matrix-store interfaces for Hub nine packs (**7B**);
- flat CpG-to-gene Deep Set and hierarchical residual-path baselines (v0.1 freezes);
- **RBS → gene cascade + direct leftover** topology (**7F**; no TBS);
- gene-only architecture selection on `explicit_only` (**7G′ Stage A** — done);
- fold-selected panel + full-model Stage B (**7G′** — GPU pending);
- static CpGPT sequence-adapter features exported offline;
- optional MethylGPT CpG-token priors as an ablation;
- masked age / tissue / sex (and Hub disease/cancer hygiene in **7E′**) heads;
- study-grouped OOF cross-fitting after 7G′ (Milestone **7**);
- array missingness and manifest-downsampling tests.

Stage 0 deliberately excludes dynamic foundation-model token extraction, LoRA,
imputation during training, episignature classifiers, epivariant calling,
production long-read training, ClickHouse, and default TileDB migration.

## Server layout

All durable and transient files must remain under `/data`. Defaults are
project-local so bootstrap does not need shared `/data/datasets` ownership:

```text
$MBS_ROOT                 Git working tree
$MBS_DATA_ROOT            data/ (canonical + staging)
$MBS_SCRATCH_ROOT         scratch/ (temporary computation)
$MBS_CACHE_ROOT           cache/ (project + tool caches)
$MBS_ARTIFACT_ROOT        artifacts/ (runs, checkpoints, scores)
```

Do not place datasets, environments, checkpoints, model weights, or caches under
`$HOME`.

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

The repository does not contain research data, pretrained weights, or copies of
the reference repositories. See [`docs/WORKSPACE.md`](docs/WORKSPACE.md) and
[`scripts/add_reference_submodules.sh`](scripts/add_reference_submodules.sh).

## Architecture

Product score path (after **7F**):

```text
Observed CpGs
  -> typed regions → RBS (one score per region)
       -> gene-allocated RBS → pool → MBS[s, g]
  -> orphan multi-CpG regions → orphan RBS columns
  -> remaining CpGs → direct (not tiled)
```

Locked Stage A gene encoder (**`P2-G`**): CascadeDeepSet on gene-linked CpGs
only (`gene_allocation: explicit_only`), **max** CpG→region and **max**
region→gene, 15 epochs. Primary metric: test-only **`mbs_e2e`**.

DeepRVAT-compatible flat baseline (still retained):

```text
CpG -> shared phi -> elementwise max by gene -> shared rho -> sigmoid MBS
```

Phenotype heads operate on the gene-score vector; exported association MBS
follows the orientation contract ([ADR 0008](docs/adr/0008-score-identifiability.md)).

End-to-end docs: [`docs/SCORING_PIPELINE.md`](docs/SCORING_PIPELINE.md),
[`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md),
[`docs/ANNOTATION_GRAPH.md`](docs/ANNOTATION_GRAPH.md).
Probe coverage: [`docs/PROBE_ANNOTATION_COVERAGE.md`](docs/PROBE_ANNOTATION_COVERAGE.md).
Data inventory: [`docs/DATA_CATALOG.md`](docs/DATA_CATALOG.md).

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
uv run mbs catalog refresh-release
uv run mbs catalog validate-release
uv run mbs inspect ewas-metadata
uv run mbs train flat --overfit-fixture

# 7G′ Stage A gene-only screen (GPU host)
# CUDA_VISIBLE_DEVICES=0 uv run python scripts/run_7g_gene_only_probe.py \
#   --config configs/experiment/stage0_7g_gene_only_probe.yaml --device cuda

# 7G′ Stage B (after Stage A lock)
# CUDA_VISIBLE_DEVICES=0 uv run python scripts/run_7g_prime_stage_b.py --device cuda
```

### Live monitoring (TensorBoard + TUI)

When `logging.tensorboard: true`, **`mbs train flat` starts TensorBoard by
default**. In a second SSH session:

```bash
source scripts/activate_data_environment.sh
uv run mbs monitor --run-id <run-id>
```

Browser over SSH: `ssh -L 6006:localhost:6006 <user>@<host>` →
http://localhost:6006. Details:
[`docs/EXPERIMENT_PROTOCOL.md`](docs/EXPERIMENT_PROTOCOL.md).

Hub/Atlas metadata: [`docs/EWAS_METADATA.md`](docs/EWAS_METADATA.md).  
Hub downloads: [`docs/EWAS_DATA.md`](docs/EWAS_DATA.md).

## Repository policy

Committed:

- Python source;
- SQL schemas and views;
- YAML configurations;
- documentation and architecture decisions;
- small synthetic test fixtures;
- artifact manifests, checksums, and inspection reports (not raw matrices).

Never committed:

- sample-level methylation data;
- IDAT, BAM, VCF, Arrow, Parquet, Zarr, HDF5, or SQLite data artifacts;
- pretrained checkpoints or embeddings;
- secrets and credentials;
- generated run outputs under `artifacts/` (except documented small reports);
- Docker layer data.

## Status

Milestones **1–7G′ Stage A** are done. Authoritative checklist:
[`docs/TODO_PIPELINE.md`](docs/TODO_PIPELINE.md).

| Done | What shipped |
|------|----------------|
| Annotation + static features | GRCh38 graphs (five-role + CGI/tile v2); offline CpGPT locus features |
| Pilot + Hub matrices | GSE35069; nine Hub full packs; ATS GSM-union 13 548 |
| Flat / hier deepMAT v0.1 | Frozen phenotype baselines (do not overwrite) |
| 7A–7E′ | Release + census; architecture corrections; Level-1 MAD; 3×2 CV; Hub multitask hygiene |
| **7F** | RBS→gene cascade + direct leftover; **no TBS** ([ADR 0009](docs/adr/0009-drop-tbs-scores.md)) |
| **7G** | Methylation-only full eval; classical vs cascade on ATS folds |
| **7G′ Stage A** | Gene-only `explicit_only` panel; test-only `mbs_e2e`; DeepRVAT screen (pooling / vector / one-hop / annotations) |

**Trustworthy Stage A numbers** (`explicit_only`, test split):

| Arm | Tissue macro-F1 | Notes |
|-----|----------------:|-------|
| `C-mvalue-enet-G` | **0.388** | Classical leader on same 51 375 CpGs |
| `P2-G` `mbs_enet` | 0.385 | Frozen MBS + elastic-net heads |
| `P2-G` `mbs_e2e` | **0.373** | **Locked cascade** (max/max, 15 ep) |
| Screen / one-hop | ≤0.359 / ~0.12 | No Tier-2 promote; prefer M-only annotations |

Report:
[`reports/inspection/stage0_7g_gene_only_probe/analysis.md`](reports/inspection/stage0_7g_gene_only_probe/analysis.md).

**Next — 7G′ Stage B GPU:** fold-safe `C-mvalue-enetS`, `N-cascade-S`,
`N-light-type`, post-hoc fusion arms, `direct_cpg.zarr`. Runner:
`scripts/run_7g_prime_stage_b.py`. Milestone **7** OOF starts only after Stage B.

Public model name remains **deepMAT**; package/CLI stay `mbs` /
`methyl-burden-score`.

## Licensing

A project source-code license has not yet been selected. Reference papers,
source repositories, datasets, and pretrained model weights have separate
licenses and must be reviewed independently before redistribution or production
use.
