# Project map

## Start here

```text
README.md                    project scope and quick start
AGENTS.md                    authoritative coding-agent rules
docs/TODO_PIPELINE.md        Stage 0 scientific milestone checklist (agents update)
docs/STRATEGIC_PLAN.md       long-term data + multimodal vision (post–Stage 0)
docs/WORKSPACE.md            server and /data layout
docs/adr/                    architecture decision records (incl. 0002 EWAS Hub)
docs/ARCHITECTURE.md         model contracts
docs/DATA_CONTRACT.md        canonical data contracts
docs/ANNOTATION_GRAPH.md     biological topology
docs/plans/                  milestone build plans (e.g. Milestone 2 graph)
docs/STATIC_FEATURES.md      CpGPT and MethylGPT artifacts
docs/EXPERIMENT_PROTOCOL.md  evaluation and controls
docs/DATA_INSPECTION.md      source acceptance workflow
docs/EWAS_DATA.md            primary open Data Hub + Atlas downloads
```

## Source package

```text
src/mbs/__init__.py      package version
src/mbs/cli.py           `mbs` command-line application
src/mbs/paths.py         /data-only filesystem policy
src/mbs/catalog.py       DuckDB catalog builder / init
src/mbs/inspect_source.py shallow source inventory reports
src/mbs/inspect_cpgcorpus.py CpGCorpus GSE/GPL scientific inspection
src/mbs/batch.py         ragged batch contract
src/mbs/segment_ops.py   permutation-invariant segment reductions
src/mbs/models.py        flat/hierarchical scorers and linear heads
src/mbs/annotation/      Stage 0 locus registry + five-role graph builder
```

These files form the current executable Stage 0 scaffold. Production ingestion,
feature stores, samplers, cross-fitting, and training orchestration are
intentionally not represented as finished modules yet.

## SQL and schemas

```text
sql/001_schema.sql                       catalog tables
sql/010_views.sql                        inspection views
schemas/matrix_manifest.schema.json      canonical matrix manifest
schemas/graph_manifest.schema.json       annotation graph manifest
schemas/static_feature_manifest.schema.json static feature manifest
```

The SQL and JSON schemas are normative interfaces. Python code should validate against them rather than creating undocumented columns.

## Configurations

```text
configs/experiment/stage0_flat.yaml      exact DeepRVAT-style baseline
configs/experiment/stage0_hier_max.yaml  hierarchical reference model
configs/local/                           machine-specific overrides, ignored
```

A training run copies its fully resolved configuration into the run artifact directory.

## Tests

```text
tests/unit/test_paths.py           filesystem policy
tests/unit/test_catalog_schema.py  real SQL schema init
tests/unit/test_cli.py             doctor / catalog / inspect CLI
tests/unit/test_segment_ops.py     reductions and permutation invariance
tests/unit/test_models.py          model and missing-gene invariants
tests/integration/test_smoke.py    catalog plus model smoke test
```

Tiny synthetic fixtures belong under `tests/fixtures/`. Never depend on private or large data in CI.

## Reference sources

```text
vendor/SOURCES.lock.yaml     reviewed repositories and commits
vendor/README.md             submodule policy
vendor/epicv2_manifest       EPICv2 reannotation code (Zenodo -> data/raw/manifests)
vendor/infinium_annotation   Zhou-lab Infinium probe coords/masks (Milestone 2)
vendor/methylcapsnet         capsule/region taxonomy reference (Milestone 2)
vendor/*                     read-only submodule pointers after installation
docs/EWAS_DATA.md            EWAS Atlas + DataHub download inventory (primary)
docs/CPGCORPUS_STAGE0.md     optional/historical CpGCorpus + labeling GSE list
reports/inspection/raw_inventory/  sizes, schemas, top rows for organizing ingest
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

## Stage 0 CLI surface

```text
mbs doctor                 validate /data paths and environment
mbs catalog init           create dirs + apply sql/*.sql schema
mbs catalog build          rebuild catalog with explicit paths
mbs inspect source         shallow raw-source inventory report
mbs inspect cpgcorpus-gpl  GSE/GPL layout, alignment, beta QC, metadata
```

## Planned next modules

Recommended implementation order:

```text
src/mbs/ingest/
src/mbs/annotations/
src/mbs/stores/
src/mbs/data/
src/mbs/training/
src/mbs/evaluation/
tools/cpgpt_export/
tools/methylgpt_export/
```

Do not create all modules as empty placeholders. Add one when its contract and tests are ready. Deeper inspection helpers may grow beside `inspect_source.py` when needed.

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