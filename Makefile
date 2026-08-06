SHELL := /usr/bin/env bash
.ONESHELL:
.SHELLFLAGS := -euo pipefail -c

PROJECT_ROOT ?= $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
DATA_ROOT ?= $(PROJECT_ROOT)/data
SCRATCH_ROOT ?= $(PROJECT_ROOT)/scratch
CACHE_ROOT ?= $(PROJECT_ROOT)/cache
ARTIFACT_ROOT ?= $(PROJECT_ROOT)/artifacts

.PHONY: help bootstrap activate doctor sync lint format typecheck test-fast test test-cov clean catalog-init catalog-build agent-context references download-cpgcorpus download-cpgcorpus-gse download-ewas-atlas download-ewas-datahub download-ewas-study download-ewas-family download-manifests download-gencode download-cpg-islands setup-methylgpt download-methylgpt export-cpgpt-static export-ewas-sample-info

help:
	@printf '%s\n' \
	  'bootstrap      Create project-local data dirs and install uv environment' \
	  'doctor         Validate environment, paths, imports, and disk space' \
	  'sync           Synchronize Python dependencies with uv' \
	  'lint           Run Ruff lint checks' \
	  'format         Format source and tests with Ruff' \
	  'typecheck      Run Pyright' \
	  'test-fast      Run unit tests' \
	  'test           Run unit and non-slow integration tests' \
	  'test-cov       Run tests with coverage' \
	  'catalog-init   Create dirs and apply sql/*.sql to the default DuckDB catalog' \
	  'catalog-build  Build the DuckDB catalog from SQL and Parquet inputs' \
	  'agent-context  Print a concise context summary for coding agents' \
	  'references     Add pinned reference repositories as submodules' \
	  'download-cpgcorpus     Sync full CpGCorpus into data/raw/cpgcorpus' \
	  'download-cpgcorpus-gse Sync Stage 0 GSE list from CpGCorpus' \
	  'download-ewas-atlas    Download EWAS Atlas batch TSVs' \
	  'download-ewas-datahub  Download EWAS DataHub methylation packs' \
	  'download-ewas-study    Download one EWAS_db study (STUDY=GSE35069)' \
	  'download-ewas-family   Download one phenotype family (FAMILY=age)' \
	  'export-ewas-sample-info Export sample-info zip → Parquet (FAMILY=tissue)' \
	  'download-manifests     Download EPICv2 Zenodo reannotated manifest' \
	  'download-gencode       Download GENCODE v38 annotation GTF' \
	  'download-cpg-islands   Download UCSC hg38 CpG island table' \
	  'setup-methylgpt        Create .venv-methylgpt and install MethylGPT deps' \
	  'download-methylgpt     Download MethylGPT medium weights + probe IDs' \
	  'export-cpgpt-static    Export CpGPT2M sequence-adapter static features' \
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

catalog-init:
	uv run mbs catalog init

catalog-build:
	uv run mbs catalog build \
	  --database "$(DATA_ROOT)/canonical/catalog/catalog.duckdb" \
	  --sql-dir "$(PROJECT_ROOT)/sql" \
	  --parquet-root "$(DATA_ROOT)/canonical/catalog/tables"

agent-context:
	bash scripts/agent_context.sh

references:
	bash scripts/add_reference_submodules.sh

download-cpgcorpus:
	bash scripts/download_cpgcorpus.sh

download-cpgcorpus-gse:
	bash scripts/download_cpgcorpus_gse.sh

download-ewas-atlas:
	bash scripts/download_ewas_atlas.sh

download-ewas-datahub:
	bash scripts/download_ewas_datahub.sh

STUDY ?= GSE35069
download-ewas-study:
	bash scripts/download_ewas_datahub_study.sh $(STUDY)

FAMILY ?= age
download-ewas-family:
	bash scripts/download_ewas_phenotype_family.sh $(FAMILY)

export-ewas-sample-info:
	uv run python -c "from pathlib import Path; from mbs.paths import DataPaths; from mbs.registry.sample_info import export_family_from_data_root; p=DataPaths.from_environment(); print(export_family_from_data_root(p.data_root, '$(FAMILY)'))"

download-manifests:
	bash scripts/download_manifests.sh

download-gencode:
	bash scripts/download_gencode.sh

download-cpg-islands:
	bash scripts/download_cpg_islands.sh

setup-methylgpt:
	bash scripts/setup_methylgpt_env.sh

download-methylgpt:
	bash scripts/download_methylgpt_weights.sh --medium

export-cpgpt-static:
	uv sync --all-groups --extra cpgpt
	uv run --extra cpgpt mbs features export-cpgpt --feature-set-id cpgpt2m_adapter_128_v1

clean:
	find src tests -type d -name __pycache__ -prune -exec rm -rf {} +
	rm -rf .pytest_cache .ruff_cache .coverage htmlcov
