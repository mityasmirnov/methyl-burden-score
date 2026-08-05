# Project map

## Start here

```text
README.md                    project scope and quick start
AGENTS.md                    authoritative coding-agent rules
docs/WORKSPACE.md            server and /data layout
docs/ARCHITECTURE.md         model contracts
docs/DATA_CONTRACT.md        canonical data contracts
docs/ANNOTATION_GRAPH.md     biological topology
docs/STATIC_FEATURES.md      CpGPT and MethylGPT artifacts
docs/EXPERIMENT_PROTOCOL.md  evaluation and controls
docs/DATA_INSPECTION.md      source acceptance workflow
```

## Source package

```text
src/mbs/__init__.py      package version
src/mbs/cli.py           `mbs` command-line application
src/mbs/paths.py         /data-only filesystem policy
src/mbs/catalog.py       DuckDB catalog builder
src/mbs/batch.py         ragged batch contract
src/mbs/segment_ops.py   permutation-invariant segment reductions
src/mbs/models.py        flat/hierarchical scorers and linear heads
```

These files form the current executable Stage 0 scaffold. Production ingestion, graph construction, feature stores, samplers, cross-fitting, and training orchestration are intentionally not represented as finished modules yet.

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
tests/unit/test_paths.py         filesystem policy
tests/unit/test_segment_ops.py   reductions and permutation invariance
tests/unit/test_models.py        model and missing-gene invariants
tests/integration/test_smoke.py  catalog plus model smoke test
```

Tiny synthetic fixtures belong under `tests/fixtures/`. Never depend on private or large data in CI.

## Reference sources

```text
vendor/SOURCES.lock.yaml  reviewed repositories and commits
vendor/README.md          submodule policy
vendor/*                  read-only submodules after local installation
```

Reference code must not become a runtime import path.

## Server scripts

```text
scripts/bootstrap_server.sh          create /data workspace and uv environment
scripts/activate_data_environment.sh export /data-only paths and caches
scripts/check_no_home_paths.sh       detect policy violations
scripts/add_reference_submodules.sh  install and commit reference submodules
scripts/agent_context.sh             concise agent context report
```

## Containers

```text
containers/Dockerfile            development image
containers/daemon.json.example   Docker daemon /data root example
compose.yaml                     bind-mounted development service
.dockerignore                    excludes data and binary artifacts
```

The host Docker daemon must already use `/data/docker` before building large images.

## Planned next modules

Recommended implementation order:

```text
src/mbs/inspect/
src/mbs/ingest/
src/mbs/annotations/
src/mbs/stores/
src/mbs/data/
src/mbs/training/
src/mbs/evaluation/
tools/cpgpt_export/
tools/methylgpt_export/
```

Do not create all modules as empty placeholders. Add one when its contract and tests are ready.

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

When a task requires raw data, provide a catalog query or sanitized inspection report instead of asking an agent to recursively inspect `/data/datasets`.