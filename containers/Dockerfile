FROM python:3.11-slim-bookworm

ARG DEBIAN_FRONTEND=noninteractive

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    MBS_ROOT=/data/projects/methyl-burden-score \
    MBS_PROJECT_ROOT=/data/projects/methyl-burden-score \
    MBS_DATA_ROOT=/data/projects/methyl-burden-score/data \
    MBS_SCRATCH_ROOT=/data/projects/methyl-burden-score/scratch \
    MBS_CACHE_ROOT=/data/projects/methyl-burden-score/cache \
    MBS_ARTIFACT_ROOT=/data/projects/methyl-burden-score/artifacts \
    MBS_DOCKER_ROOT=/data/projects/methyl-burden-score/docker \
    UV_CACHE_DIR=/data/projects/methyl-burden-score/cache/uv \
    PIP_CACHE_DIR=/data/projects/methyl-burden-score/cache/pip \
    HF_HOME=/data/projects/methyl-burden-score/cache/huggingface \
    TORCH_HOME=/data/projects/methyl-burden-score/cache/torch \
    XDG_CACHE_HOME=/data/projects/methyl-burden-score/cache/xdg \
    TMPDIR=/data/projects/methyl-burden-score/scratch/tmp

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
        "$MBS_DATA_ROOT" \
        "$MBS_SCRATCH_ROOT/tmp" \
        "$MBS_CACHE_ROOT/uv" \
        "$MBS_CACHE_ROOT/pip" \
        "$MBS_CACHE_ROOT/huggingface" \
        "$MBS_CACHE_ROOT/torch" \
        "$MBS_CACHE_ROOT/xdg" \
        "$MBS_ARTIFACT_ROOT" \
    && uv sync --all-groups --all-extras

ENV PATH="/data/projects/methyl-burden-score/.venv/bin:${PATH}"

CMD ["uv", "run", "mbs", "doctor"]
