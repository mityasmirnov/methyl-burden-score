FROM python:3.11-slim-bookworm

ARG DEBIAN_FRONTEND=noninteractive

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    MBS_PROJECT_ROOT=/data/projects/methyl-burden-score \
    MBS_DATA_ROOT=/data/datasets/methyl-burden-score \
    MBS_SCRATCH_ROOT=/data/scratch/methyl-burden-score \
    MBS_CACHE_ROOT=/data/cache/methyl-burden-score \
    MBS_ARTIFACT_ROOT=/data/artifacts/methyl-burden-score \
    MBS_DOCKER_ROOT=/data/docker \
    UV_CACHE_DIR=/data/cache/uv \
    PIP_CACHE_DIR=/data/cache/pip \
    HF_HOME=/data/cache/huggingface \
    TORCH_HOME=/data/cache/torch \
    XDG_CACHE_HOME=/data/cache/xdg \
    TMPDIR=/data/scratch/methyl-burden-score/tmp

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        ca-certificates \
        curl \
        git \
        build-essential \
    && rm -rf /var/lib/apt/lists/* \
    && python -m pip install --no-cache-dir uv

WORKDIR /data/projects/methyl-burden-score

COPY pyproject.toml README.md AGENTS.md ./
COPY src ./src
COPY tests ./tests
COPY configs ./configs
COPY schemas ./schemas
COPY sql ./sql

RUN mkdir -p \
        /data/datasets/methyl-burden-score \
        /data/scratch/methyl-burden-score/tmp \
        /data/cache/methyl-burden-score \
        /data/cache/uv \
        /data/cache/pip \
        /data/cache/huggingface \
        /data/cache/torch \
        /data/cache/xdg \
        /data/artifacts/methyl-burden-score \
    && uv sync --all-groups --all-extras

ENV PATH="/data/projects/methyl-burden-score/.venv/bin:${PATH}"

CMD ["uv", "run", "mbs", "doctor"]
