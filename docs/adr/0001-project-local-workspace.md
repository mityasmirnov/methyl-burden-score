# ADR 0001: Project-local workspace roots

## Status

Accepted

## Context

Stage 0 needs a reproducible `/data`-only filesystem policy on shared servers. An earlier draft layout used shared trees such as `/data/datasets/methyl-burden-score` and `/data/tools/uv`. That layout requires cross-project ownership and diverges from what `.env.example`, `DataPaths`, and the bootstrap scripts actually create.

## Decision

Default all MBS writable roots under the Git working tree:

```text
MBS_ROOT=/data/projects/methyl-burden-score
MBS_DATA_ROOT=$MBS_ROOT/data
MBS_SCRATCH_ROOT=$MBS_ROOT/scratch
MBS_CACHE_ROOT=$MBS_ROOT/cache
MBS_ARTIFACT_ROOT=$MBS_ROOT/artifacts
MBS_DOCKER_ROOT=$MBS_ROOT/docker
```

Install `uv` under `$MBS_ROOT/.tools/uv`. Nest tool caches (`huggingface`, `uv`, `pip`, `torch`, …) under `$MBS_CACHE_ROOT`.

The path policy still requires every configured root to be an absolute path under `/data`. Shared mounts may override defaults via `.env`, but project-local paths are the documented and tested default.

Reference repositories under `vendor/` remain Git submodule **gitlinks** (mode `160000`). Initialized submodule working trees on disk are expected after `git submodule update --init`. Agents and contributors must never copy upstream trees into the repository as ordinary tracked directories.

## Consequences

- Bootstrap works without shared `/data/datasets` group ownership.
- `docs/WORKSPACE.md`, `.env.example`, shell scripts, and `src/mbs/paths.py` stay aligned.
- Docker **daemon** root (`/data/docker`) remains a host/admin concern, separate from `$MBS_DOCKER_ROOT`.
- Sanitized inspection reports stay under repo `reports/inspection/`; bulky run outputs stay under `$MBS_ARTIFACT_ROOT`.
