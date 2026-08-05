SHELL := /usr/bin/env bash
.ONESHELL:
.SHELLFLAGS := -euo pipefail -c

PROJECT_ROOT ?= /data/projects/methyl-burden-score
DATA_ROOT ?= /data/datasets/methyl-burden-score
SCRATCH_ROOT ?= /data/scratch/methyl-burden-score
CACHE_ROOT ?= /data/cache/methyl-burden-score
ARTIFACT_ROOT ?= /data/artifacts/methyl-burden-score

.PHONY: help bootstrap activate doctor sync lint format typecheck test-fast test test-cov clean catalog-build agent-context references

help:
	@printf '%s\n' \
	  'bootstrap      Create /data directories and install uv environment' \
	  'doctor         Validate environment, paths, imports, and disk space' \
	  'sync           Synchronize Python dependencies with uv' \
	  'lint           Run Ruff lint checks' \
	  'format         Format source and tests with Ruff' \
	  'typecheck      Run Pyright' \
	  'test-fast      Run unit tests' \
	  'test           Run unit and non-slow integration tests' \
	  'test-cov       Run tests with coverage' \
	  'catalog-build  Build the DuckDB catalog from SQL and Parquet inputs' \
	  'agent-context  Print a concise context summary for coding agents' \
	  'references     Add pinned reference repositories as submodules' \
	  'clean          Remove generated local Python caches only'

bootstrap:
	bash scripts/bootstrap_server.sh

activate:
	@printf 'Run: source scripts/activate_data_environment.sh\n'

doctor:
	uv run mbs doctor

sync:
	uv sync --all-groups --all-extras

lint:
	uv run ruff check .

format:
	uv run ruff format src tests

typecheck:
	uv run pyright

test-fast:
	uv run pytest tests/unit

test:
	uv run pytest tests/unit tests/integration -m 'not slow'

test-cov:
	uv run pytest --cov=mbs --cov-report=term-missing tests/unit tests/integration -m 'not slow'

catalog-build:
	uv run mbs catalog build \
	  --database "$(DATA_ROOT)/canonical/catalog/catalog.duckdb" \
	  --sql-dir "$(PROJECT_ROOT)/sql" \
	  --parquet-root "$(DATA_ROOT)/canonical/catalog/tables"

agent-context:
	bash scripts/agent_context.sh

references:
	bash scripts/add_reference_submodules.sh

clean:
	find src tests -type d -name __pycache__ -prune -exec rm -rf {} +
	rm -rf .pytest_cache .ruff_cache .coverage htmlcov
