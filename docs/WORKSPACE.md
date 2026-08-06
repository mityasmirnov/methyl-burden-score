# Server workspace

## Fixed project location

The Git working tree lives at:

```text
/data/projects/methyl-burden-score
```

Open this exact folder in Cursor or another coding agent. Do not open `/data`, because doing so will index unrelated projects and potentially very large datasets.

## Directory layout

Defaults are **project-local** under the Git working tree. Any absolute path under `/data` still passes the path-policy check, but overrides are optional — bootstrap does not require shared `/data/datasets` ownership.

```text
/data/projects/methyl-burden-score/          # $MBS_ROOT
├── data/                                    # $MBS_DATA_ROOT
│   ├── raw/                                 # immutable source downloads
│   │   ├── cpgcorpus/
│   │   ├── ewas_atlas/
│   │   ├── ewas_datahub/
│   │   ├── manifests/
│   │   └── methylgpt/           # pretrained checkpoints + probe vocab
│   ├── staging/                             # disposable conversions
│   └── canonical/
│       ├── catalog/
│       │   ├── catalog.duckdb
│       │   └── tables/                      # Parquet metadata tables
│       ├── matrices/
│       ├── annotations/
│       ├── graphs/
│       └── static_features/
├── scratch/                                 # $MBS_SCRATCH_ROOT (TMPDIR)
│   ├── tmp/
│   └── downloads/
├── cache/                                   # $MBS_CACHE_ROOT
│   ├── huggingface/
│   ├── uv/
│   ├── pip/
│   ├── torch/
│   └── …
├── artifacts/                               # $MBS_ARTIFACT_ROOT
│   ├── runs/
│   ├── checkpoints/
│   ├── scores/
│   ├── reports/
│   ├── logs/downloads/
│   └── wandb/
├── docker/                                  # $MBS_DOCKER_ROOT (local bind/mount helper)
├── .tools/uv/                               # project-local uv install
├── .venv/                                   # main MBS / optional CpGPT env
└── .venv-methylgpt/                         # MethylGPT-only (torchtext-compatible)
```

Environment variables (see `.env.example`):

```text
MBS_ROOT            /data/projects/methyl-burden-score
MBS_DATA_ROOT       $MBS_ROOT/data
MBS_SCRATCH_ROOT    $MBS_ROOT/scratch
MBS_CACHE_ROOT      $MBS_ROOT/cache
MBS_ARTIFACT_ROOT   $MBS_ROOT/artifacts
MBS_DOCKER_ROOT     $MBS_ROOT/docker
```

`MBS_PROJECT_ROOT` remains a compatibility alias for `MBS_ROOT`.

The host Docker **daemon** data root is separate and typically `/data/docker` (see below). That is an administrator concern, not the project-local `$MBS_DOCKER_ROOT` helper directory.

## Initial server setup

```bash
mkdir -p /data/projects
cd /data/projects
git clone git@github.com:mityasmirnov/methyl-burden-score.git
cd methyl-burden-score

cp .env.example .env
source scripts/activate_data_environment.sh
bash scripts/bootstrap_server.sh
```

The bootstrap script installs `uv` under `$MBS_ROOT/.tools/uv/bin`, creates `.venv` inside the project, and places all tool caches under `$MBS_CACHE_ROOT`.

**PyTorch / CUDA:** On Linux and Windows, `uv sync` installs torch from the
PyTorch **cu128** wheel index (matches this server's NVIDIA driver 570.x /
CUDA 12.8). Confirm GPUs with:

```bash
uv run python -c "import torch; print(torch.__version__, torch.cuda.is_available(), torch.cuda.device_count())"
```

Expect something like `2.11.0+cu128 True 3`. If you see `+cu130` and
`cuda.is_available() == False`, re-sync from a clean lock that still points at
`https://download.pytorch.org/whl/cu128` (see `pyproject.toml`).

Optional foundation-model export tooling is **not** in the default sync.

**CpGPT** (optional extra on the main `.venv`):

```bash
uv sync --all-groups --extra cpgpt
```

Weights and human dependencies go under `$MBS_CACHE_ROOT/huggingface` via
`download_cpgpt` (see `docs/STATIC_FEATURES.md`). Never under `$HOME/.cache`.

**MethylGPT** needs a separate env (torchtext ABI vs main torch). After bootstrap:

```bash
make setup-methylgpt
make download-methylgpt
```

Checkpoints land under `$MBS_DATA_ROOT/raw/methylgpt/` (see
`docs/STATIC_FEATURES.md`). Do not install MethylGPT into `.venv`.

## Shell startup

Do not add large environment initializations to `$HOME`. A small shell line that sources a script under `/data` is acceptable when desired:

```bash
source /data/projects/methyl-burden-score/scripts/activate_data_environment.sh
```

For project-specific use, source this only after entering the project.

## Cursor and coding agents

Open:

```text
/data/projects/methyl-burden-score
```

Agent-visible material:

- `src/`, `tests/`, `configs/`, `schemas/`, and `sql/`;
- project documentation and ADRs under `docs/adr/`;
- sanitized reports in `reports/inspection/`;
- targeted files in read-only reference submodules.

Agent-excluded material:

- raw and canonical sample matrices;
- caches, checkpoints, and model weights;
- Zarr, HDF5, BAM, IDAT, SQLite, and DuckDB binaries;
- private sample metadata;
- bulk vendor submodule trees (use targeted reads only).

Run this before handing a task to an agent:

```bash
make agent-context
```

## Docker under `/data`

The repository includes `containers/daemon.json.example`. Moving the Docker daemon root requires administrative access and should be performed only after checking the current daemon configuration.

Typical administrator procedure:

```bash
sudo systemctl stop docker
sudo mkdir -p /data/docker
sudo rsync -aHAX --numeric-ids /var/lib/docker/ /data/docker/
sudo install -D -m 0644 containers/daemon.json.example /etc/docker/daemon.json
sudo systemctl start docker
sudo docker info --format '{{.DockerRootDir}}'
```

Expected output:

```text
/data/docker
```

Do not delete `/var/lib/docker` until all containers, images, volumes, and build caches have been verified under the new root. Merge the example JSON with any existing daemon settings rather than overwriting them blindly.

## Data lifecycle

### Raw

`$MBS_DATA_ROOT/raw/` is immutable. Every source file receives:

- retrieval timestamp;
- source URL or accession;
- byte size;
- SHA-256 checksum;
- data-use and license notes.

### Staging

`$MBS_DATA_ROOT/staging/` is disposable. Any artifact here must be reproducible from raw files and committed code.

### Canonical

A canonical release is immutable after it is used for a reported run. Each release contains manifests that identify its input checksums, graph version, locus registry, and conversion code commit.

### Artifacts

Every run directory under `$MBS_ARTIFACT_ROOT/runs/` contains a resolved configuration, Git commit, environment summary, split manifest, model checkpoints, metrics, and output checksums.

## Permissions

Do not make research data world-readable. A typical shared-project policy is:

```bash
umask 007
find "$MBS_DATA_ROOT" -type d -exec chmod 2770 {} +
find "$MBS_DATA_ROOT" -type f -exec chmod 0660 {} +
```

Choose the correct institutional Unix group before applying recursive permission changes.

## Space checks

```bash
df -h /data
du -sh /data/projects/methyl-burden-score
du -sh "$MBS_DATA_ROOT" "$MBS_CACHE_ROOT" "$MBS_ARTIFACT_ROOT" 2>/dev/null
docker system df
```

Run `scripts/check_no_home_paths.sh` after installing new tooling or foundation-model dependencies.
