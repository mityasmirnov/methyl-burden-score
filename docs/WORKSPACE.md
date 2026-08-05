# Server workspace

## Fixed project location

The Git working tree lives at:

```text
/data/projects/methyl-burden-score
```

Open this exact folder in Cursor or another coding agent. Do not open `/data`, because doing so will index unrelated projects and potentially very large datasets.

## Directory layout

```text
/data/
├── projects/
│   └── methyl-burden-score/          # Git repository and .venv
├── datasets/
│   └── methyl-burden-score/
│       ├── raw/                       # immutable source downloads
│       ├── staging/                   # disposable conversions
│       └── canonical/
│           ├── catalog/
│           ├── matrices/
│           ├── annotations/
│           ├── graphs/
│           └── static_features/
├── scratch/
│   └── methyl-burden-score/           # temporary jobs and TMPDIR
├── cache/
│   ├── methyl-burden-score/
│   ├── huggingface/
│   ├── uv/
│   ├── pip/
│   └── torch/
├── artifacts/
│   └── methyl-burden-score/
│       ├── runs/
│       ├── checkpoints/
│       ├── scores/
│       └── reports/
├── tools/
│   └── uv/
└── docker/                            # Docker daemon data root
```

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

The bootstrap script installs `uv` under `/data/tools/uv`, creates `.venv` inside the project, and places all caches under `/data`.

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
- project documentation;
- sanitized reports in `reports/inspection/`;
- targeted files in read-only reference submodules.

Agent-excluded material:

- raw and canonical sample matrices;
- caches, checkpoints, and model weights;
- Zarr, HDF5, BAM, IDAT, SQLite, and DuckDB binaries;
- private sample metadata.

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

`raw/` is immutable. Every source file receives:

- retrieval timestamp;
- source URL or accession;
- byte size;
- SHA-256 checksum;
- data-use and license notes.

### Staging

`staging/` is disposable. Any artifact here must be reproducible from raw files and committed code.

### Canonical

A canonical release is immutable after it is used for a reported run. Each release contains manifests that identify its input checksums, graph version, locus registry, and conversion code commit.

### Artifacts

Every run directory contains a resolved configuration, Git commit, environment summary, split manifest, model checkpoints, metrics, and output checksums.

## Permissions

Do not make research data world-readable. A typical shared-project policy is:

```bash
umask 007
find /data/datasets/methyl-burden-score -type d -exec chmod 2770 {} +
find /data/datasets/methyl-burden-score -type f -exec chmod 0660 {} +
```

Choose the correct institutional Unix group before applying recursive permission changes.

## Space checks

```bash
df -h /data
du -sh /data/projects/methyl-burden-score
du -sh /data/datasets/methyl-burden-score
du -sh /data/cache/* 2>/dev/null | sort -h
docker system df
```

Run `scripts/check_no_home_paths.sh` after installing new tooling or foundation-model dependencies.