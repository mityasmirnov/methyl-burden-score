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
14. Do not launch final 5-fold × 6-restart OOF cross-fitting until Milestones 7A–7E′ are done ([ADR 0007](docs/adr/0007-crossfit-prerequisites.md)). Do not retrain frozen v0.1 runs; 7A–7D are done and the coding gate is **7E** (graph-v2 on disk; multi-path unblocked). Hub-wide disease/cancer heads and catalog hygiene are **7E′** (unknown labels stay unknown, not controls).
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