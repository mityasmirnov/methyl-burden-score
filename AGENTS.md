# Coding-agent instructions

This file is authoritative for Cursor, Claude, Codex, and other coding agents working in this repository.

## Mission

Implement Stage 0 of the Methylation Burden Score project: a reproducible,
study-grouped, cross-fitted pipeline that maps variable sets of observed CpGs to
typed regulatory regions and then to gene-level **MBS** scores (plus, from
Milestone 7C, optional non-gene **RBS**, intergenic **TBS**, and direct CpG
contributions) for association and prediction. Current gate and ordering:
[`docs/TODO_PIPELINE.md`](docs/TODO_PIPELINE.md).

## Filesystem policy

All project files, datasets, environments, caches, temporary files, Docker data, checkpoints, and model artifacts must reside under `/data`.

Canonical locations (project-local defaults; override via `.env`):

```text
$MBS_ROOT                 Git working tree (/data/projects/methyl-burden-score)
$MBS_DATA_ROOT            data/ under the project
$MBS_SCRATCH_ROOT         scratch/ under the project
$MBS_CACHE_ROOT           cache/ under the project (incl. huggingface, uv, pip, torch)
$MBS_ARTIFACT_ROOT        artifacts/ under the project
```

Never create project artifacts in `$HOME`, `~/.cache`, `/tmp`, or `/var/lib/docker` unless an administrator has explicitly configured those paths as symlinks or mounts under `/data`.

## Editable and read-only areas

Editable:

```text
src/
tests/
configs/
sql/
docs/
schemas/
scripts/
tools/
containers/
```

Read-only references:

```text
vendor/
```

Do not modify vendored repositories from this project. Changes to a reference project must be made in its own fork and pinned here by commit SHA.

## Never inspect recursively

Do not recursively read, index, or summarize:

```text
$MBS_DATA_ROOT/
$MBS_CACHE_ROOT/
$MBS_SCRATCH_ROOT/
$MBS_ARTIFACT_ROOT/
*.zarr/
*.h5
*.hdf5
*.bam
*.cram
*.idat
*.ckpt
*.pt
*.pth
```

Use the project inspection commands, DuckDB queries, manifests, and sanitized reports under `reports/inspection/` instead.

## Scientific invariants

1. Sample IDs, study IDs, GSE IDs, donor IDs, and platform IDs must not enter the shared MBS encoder as predictive features.
2. Any normalization or phenotype-derived feature selection must be fitted inside the training fold.
3. No sample may be scored by a model trained on that sample, its donor, its technical replicate, or its held-out study group.
4. CpG and region aggregation must be permutation invariant.
5. Empty genes must produce a neutral score and a false presence mask; missingness must never be interpreted as low burden.
6. Static foundation-model features must have a manifest containing repository commit, checkpoint hash, vocabulary or locus-table hash, dimensions, dtype, genome build, and export command.
7. Graph, matrix, split, and checkpoint artifacts must be content-addressed or accompanied by checksums.
8. The exact flat DeepRVAT-style max-pooling model must remain available as a reference baseline.
9. Production logic belongs in `src/`; notebooks may only call package functions.
10. Do not advertise placeholder CLI commands as implemented.
11. Do not store sample×CpG observations as a long DuckDB table; catalog stays independent of the matrix backend ([ADR 0005](docs/adr/0005-catalog-matrix-independence.md)).
12. Missing disease/cancer labels are **unknown**, not automatic controls, unless pack documentation establishes controls.
13. Do not assign intergenic loci to the nearest gene; do not treat one-scalar residual compression as proof that noncoding CpGs are uninformative ([ADR 0006](docs/adr/0006-multipath-noncoding-scores.md)).
14. Do not launch final 5-fold × 6-restart OOF cross-fitting until Milestones
    **7F and 7G** are done ([ADR 0007](docs/adr/0007-crossfit-prerequisites.md),
    [ADR 0009](docs/adr/0009-drop-tbs-scores.md)). Do not retrain frozen v0.1
    runs; 7A–7F and **7E′** are done and the coding gate is **7G**. Hub-wide
    disease/cancer heads and catalog hygiene (**7E′**) keep unknown labels
    unknown, not controls. Product scores after 7F: gene-aggregated RBS (MBS),
    orphan RBS, and direct — **no TBS**.
15. Before averaging OOF scores, apply a score-orientation anchor; predictive MBS is not a constraint/LOEUF score ([ADR 0008](docs/adr/0008-score-identifiability.md)).

## Coding standards

- Python 3.11 or newer.
- Type annotations on public functions and dataclasses.
- Prefer standard PyTorch operations over mandatory compiled extensions.
- Use `pathlib.Path` rather than string path concatenation.
- Use LayerNorm or RMSNorm in the shared encoder; do not introduce BatchNorm without an explicit ADR.
- Keep modules small and testable.
- Raise explicit errors on manifest, checkpoint, vocabulary, coordinate, or genome-build mismatches.
- Never silently clip methylation values or drop samples during ingestion.

## Required checks

Before completing a code change, run:

```bash
uv run ruff check .
uv run ruff format --check .
uv run pyright
uv run pytest tests/unit
uv run pytest tests/integration -m "not slow"
```

For a documentation-only change, at minimum run:

```bash
uv run ruff check .
uv run pytest tests/unit
```

## Source-of-truth hierarchy

1. Architecture decision records in `docs/adr/`.
2. Schemas in `schemas/` and SQL in `sql/`.
3. Resolved experiment configuration.
4. Python implementation.
5. Notebooks and exploratory reports.

When these disagree, stop and report the inconsistency rather than guessing.

## Cursor Cloud specific instructions

This project is a pure Python CLI / library / batch pipeline (`mbs` command); there
is no web app or long-running service. TensorBoard is the only server, started
on demand during real training. **On machines with a GPU**, always pass
`--device cuda` (or `device.torch_device: cuda` in YAML) for real Hub matrix
training — see [`docs/SCORING_PIPELINE.md`](docs/SCORING_PIPELINE.md) and
[`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md). The Cursor Cloud VM has **no
GPU**, so `torch.cuda.is_available()` is `False` and everything runs on CPU
(fixtures fall back automatically; pass `--device cpu` for real configs on Cloud
only).

### /data policy on the Cloud VM (required before running anything)

`src/mbs/paths.py` requires `MBS_ROOT` and all data/cache/scratch/artifact roots
to resolve under `/data`, but the git working tree lives at `/workspace`. Setup
bind-mounts the working tree at `/data/projects/methyl-burden-score` (a bind
mount, not a symlink, so `Path.resolve()` stays under `/data`; the path-policy
checks and `tests/unit/test_cli.py` scratch fixtures depend on this). The bind
mount is runtime state and is re-established by the startup update script; if it
is ever missing, recreate it (idempotent, needs sudo):

```bash
bash scripts/cloud_data_setup.sh
```

Before running the CLI or tests, always activate from the /data project root:

```bash
cd /data/projects/methyl-burden-score
source scripts/activate_data_environment.sh   # exports MBS_* under /data + PYTHONPATH=$MBS_ROOT/src
```

Running `mbs`/`pytest` from `/workspace` without the bind mount + activate makes
path-policy checks and ~20 path/CLI tests fail with `Path policy failure` — that
is a missing-activation symptom, not a code regression.

### Dependencies

`uv` lives at `~/.local/bin` (on `PATH` via `~/.bashrc`). The startup update
script runs `uv sync --all-groups --extra training --extra analysis --frozen`.
Notes:
- Always use `--frozen`. The optional `cpgpt` extra maps to the empty vendored
  submodule `vendor/cpgpt` (`[tool.uv.sources]`), which breaks a non-frozen
  resolve. Keep `--frozen` and do **not** pass `--extra cpgpt` unless the
  submodule is populated (`git submodule update --init vendor/cpgpt`).
- The `training` extra (lightning/scikit-learn/tensorboard/torchmetrics/tqdm) is
  required even for `mbs train flat --overfit-fixture` (`loop.py` imports
  scikit-learn and `torch.utils.tensorboard`).
- After a dependency change, re-run `uv sync ... --frozen`; when the running
  `mbs` process is a training loop it does not hot-reload, restart it.

### Checks that are enforced vs. pre-existing noise

CI (`.github/workflows/ci.yml`) hard-fails only on `uv run ruff check --select
E4,E7,E9,F .`, `uv run pytest tests/unit`, and `uv run pytest tests/integration
-m "not slow"`; those pass cleanly here. The full `uv run ruff check .`, `uv run
ruff format --check .`, and `uv run pyright` (strict) report pre-existing
findings on `main` (newer ruff/pyright than the tree was written against) and CI
marks format + pyright `continue-on-error`. Do not treat those pre-existing
findings as regressions from your change.

### Real end-to-end pipeline (beyond fixtures)

`graph build`, `matrix convert*`, `features export-cpgpt`, and real training need
external downloads (`make download-*`) and populated vendored submodules
(`scripts/add_reference_submodules.sh` / `git submodule update --init`), none of
which ship in the repo. The synthetic `--overfit-fixture` paths and
`tests/fixtures/` exercise the model code without any real data.