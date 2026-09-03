# Project map

## Start here

```text
README.md                    project scope and quick start
AGENTS.md                    authoritative coding-agent rules
docs/TODO_PIPELINE.md        Stage 0 scientific milestone checklist (agents update)
docs/STRATEGIC_PLAN.md       long-term data + multimodal vision
docs/WORKSPACE.md            server and /data layout
docs/adr/                    architecture decision records (0001 workspace, 0002 EWAS Hub, 0003 Milestone 5b, 0004 unmapped retention, 0005 catalog/storage, 0006 multi-path scores, 0007 cross-fit prerequisites, 0008 score identifiability)
docs/ARCHITECTURE.md         model contracts (public name: deepMAT; package: mbs)
docs/DATA_CONTRACT.md        canonical data contracts
docs/ANNOTATION_GRAPH.md     biological topology
docs/plans/                  milestone build plans (incl. post-v0 programme)
docs/STATIC_FEATURES.md      CpGPT and MethylGPT artifacts
docs/EXPERIMENT_PROTOCOL.md  evaluation, controls, train monitoring (TB/JSONL)
docs/DATA_INSPECTION.md      source acceptance workflow
docs/EWAS_DATA.md            primary open Data Hub + Atlas downloads
docs/EWAS_METADATA.md        Atlas small tables + Hub sample-info contracts
docs/DATA_CATALOG.md         Hub packs, Ns, freeze tags
docs/SCORING_PIPELINE.md     CpG → MBS / heads narrative
```

## Source package

```text
src/mbs/__init__.py      package version
src/mbs/cli.py           `mbs` command-line application
src/mbs/paths.py         /data-only filesystem policy
src/mbs/catalog.py       DuckDB catalog builder / init (schema; 7A populates)
src/mbs/inspect_source.py shallow source inventory reports
src/mbs/inspect_cpgcorpus.py CpGCorpus GSE/GPL scientific inspection
src/mbs/inspect_ewas_metadata.py Atlas small tables + Hub sample-info profiles
src/mbs/geo_metadata.py      GEO family SOFT parse + EWAS_db-only sample backfill
src/mbs/registry/sample_info.py Hub sample-info txt/zip → Parquet export
src/mbs/batch.py         ragged batch contract
src/mbs/segment_ops.py   permutation-invariant segment reductions
src/mbs/models.py        flat/hierarchical scorers and linear heads
src/mbs/annotation/      Stage 0 locus registry + five-role / graph-v2 builder
src/mbs/static_features/ offline CpGPT sequence-adapter export + artifact I/O
src/mbs/matrix/          canonical matrix conversion (EWAS_db + Hub packs)
src/mbs/scoring/         score orientation / score_manifest (ADR 0008)
src/mbs/training/        flat / hierarchical / branch / controls / sampler
src/mbs/registry/        phenotype / source dataset registry (Milestone 5b)
src/mbs/evaluation/      metrics + study-grouped split helpers (Milestone 5b)
```

Milestones **1–6**, **7A**, and **7C** (fixture acceptance) are done.
**Current gate: 7B** (convert remaining Hub packs to canonical matrices +
inspection report). **7C** trainer/graph-v2/direct/branch fixtures are closed;
residual wiring listed in
[`plans/milestone-7c-supervised-architecture.md`](plans/milestone-7c-supervised-architecture.md).
Final OOF cross-fitting is Milestone **7**, blocked until 7A–7E — see
[`TODO_PIPELINE.md`](TODO_PIPELINE.md),
[ADR 0007](adr/0007-crossfit-prerequisites.md),
[ADR 0008](adr/0008-score-identifiability.md), and
[`plans/post-v0-scientific-programme.md`](plans/post-v0-scientific-programme.md).
**Do not retrain v0.1**.

## SQL and schemas

```text
sql/001_schema.sql                       catalog tables
sql/002_provenance_lanes.sql             provenance lanes
sql/010_views.sql                        inspection views
sql/011_census_views.sql                 phenotype census / eligibility views (7A)
schemas/release_manifest.schema.json     versioned data release manifest (7A)
schemas/matrix_manifest.schema.json      canonical matrix manifest
schemas/graph_manifest.schema.json       annotation graph manifest
schemas/static_feature_manifest.schema.json static feature manifest
schemas/phenotype_registry.schema.json   phenotype / source registry (5b)
schemas/sample_phenotype_table.schema.json unified sample×task table (5c)
```

The SQL and JSON schemas are normative interfaces. Python code should validate against them rather than creating undocumented columns.

## Configurations

```text
configs/experiment/stage0_flat.yaml        exact DeepRVAT-style baseline
configs/experiment/stage0_flat_pilot.yaml  GSE35069 cell-type pilot train config
configs/experiment/stage0_flat_age_holdout.yaml   age study-holdout fixture config
configs/experiment/stage0_flat_tissue_holdout.yaml tissue study-holdout fixture
configs/experiment/stage0_flat_multitask.yaml     Hub multitask (Milestone 5c)
configs/experiment/stage0_flat_deeprvat_full.yaml max-N age/tissue/sex (5d)
configs/experiment/stage0_hier_deeprvat_full.yaml hierarchical v0.1 (6)
configs/experiment/stage0_hier_max.yaml    hierarchical sketch config
configs/data/phenotype_registry.yaml     dataset registry (Milestone 5b)
configs/local/                             machine-specific overrides, ignored
```

A training run copies its fully resolved configuration into the run artifact directory.

## Tests

```text
tests/unit/test_paths.py           filesystem policy
tests/unit/test_catalog_schema.py  real SQL schema init
tests/unit/test_cli.py             doctor / catalog / inspect CLI
tests/unit/test_segment_ops.py     reductions and permutation invariance
tests/unit/test_models.py          model and missing-gene invariants
tests/unit/test_training_flat.py   flat baseline phenotypes / overfit / CLI
tests/unit/test_phenotype_registry.py  registry loader / checksums
tests/unit/test_evaluation.py      metrics + study-grouped splits
tests/integration/test_smoke.py    catalog plus model smoke test
```

Tiny synthetic fixtures belong under `tests/fixtures/`. Never depend on private or large data in CI.

## Reference sources

```text
vendor/SOURCES.lock.yaml     reviewed repositories and commits
vendor/README.md             submodule policy
vendor/epicv2_manifest       EPICv2 reannotation code (Zenodo -> data/raw/manifests)
vendor/infinium_annotation   Zhou-lab Infinium probe coords/masks (Milestone 2)
vendor/methylcapsnet         capsule/region taxonomy reference (Milestone 2 / 7C)
vendor/*                     read-only submodule pointers after installation
docs/EWAS_DATA.md            EWAS Atlas + DataHub download inventory (primary)
docs/EWAS_METADATA.md        Atlas small tables + Hub sample-info contracts
docs/CPGCORPUS_STAGE0.md     optional/historical CpGCorpus + labeling GSE list
reports/inspection/raw_inventory/  sizes, schemas, top rows for organizing ingest
reports/inspection/ewas_metadata_structure/  Atlas + sample-info structure profiles
reports/inspection/ewas_datahub_samples/     unpacked Hub sample_*.txt (Cursor-visible)
```

Reference code must not become a runtime import path. Downloaded manifests and corpora belong under `$MBS_DATA_ROOT/raw/`.

## Server scripts

```text
scripts/bootstrap_server.sh              create /data workspace and uv environment
scripts/activate_data_environment.sh     export /data-only paths and caches
scripts/check_no_home_paths.sh           detect policy violations
scripts/add_reference_submodules.sh      install and commit reference submodules
scripts/agent_context.sh                 concise agent context report
scripts/download_cpgcorpus.sh            full CpGCorpus sync (explicit invoke)
scripts/download_cpgcorpus_gse.sh        Stage 0 GSE sync (explicit invoke)
scripts/download_cpgcorpus_background.sh nohup wrapper for CpGCorpus syncs
```

## Containers

```text
containers/Dockerfile            development image
containers/daemon.json.example   Docker daemon /data root example
compose.yaml                     bind-mounted development service
.dockerignore                    excludes data and binary artifacts
```

The host Docker daemon must already use `/data/docker` before building large images.

## Stage 0 CLI surface (implemented)

```text
mbs doctor                 validate /data paths and environment
mbs catalog init           create dirs + apply sql/*.sql schema
mbs catalog build          rebuild catalog with explicit paths
mbs catalog refresh-release  populate deepmat-data-v1 (+ EWAS_db upsert)
mbs catalog validate-release validate release manifest + catalog
mbs catalog phenotype-census write phenotype census report
mbs catalog trait-eligibility write trait eligibility report
mbs inspect source         shallow raw-source inventory report
mbs inspect cpgcorpus-gpl  GSE/GPL layout, alignment, beta QC, metadata
mbs inspect ewas-metadata  Atlas small tables + Hub sample-info structure
mbs graph build            locus registry + five-role annotation graph
mbs features export-cpgpt  offline CpGPT2M sequence-adapter static features
mbs matrix convert         EWAS_db study → canonical matrix store
mbs matrix convert-pack    Hub phenotype pack → canonical matrix
mbs matrix index-hub-packs Virtual multi-store index + GSM overlap check
mbs train flat             flat DeepRVAT-style baseline
mbs train hierarchical     hierarchical + residual-slot v0.1
mbs monitor                live Rich TUI for a run (metrics.jsonl + GPU + ckpts)
mbs phenotypes build-multitask-table   Hub age+tissue(+sex) phenotype table
```

Milestone **7B+** work follows [`TODO_PIPELINE.md`](TODO_PIPELINE.md)
(7A→7E→7). Do not create empty placeholders. Add one module when its contract and
tests are ready.

## Coding-agent task boundary

Before editing:

```bash
source scripts/activate_data_environment.sh
make agent-context
```

Before completion:

```bash
make lint
make typecheck
make test
```

When a task requires raw data, provide a catalog query or sanitized inspection report instead of asking an agent to recursively inspect `$MBS_DATA_ROOT`.
