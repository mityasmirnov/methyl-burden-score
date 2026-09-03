SHELL := /bin/bash
.ONESHELL:
.SHELLFLAGS := -euo pipefail -c

PROJECT_ROOT ?= $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
DATA_ROOT ?= $(PROJECT_ROOT)/data
SCRATCH_ROOT ?= $(PROJECT_ROOT)/scratch
CACHE_ROOT ?= $(PROJECT_ROOT)/cache
ARTIFACT_ROOT ?= $(PROJECT_ROOT)/artifacts

.PHONY: help bootstrap activate doctor sync lint format typecheck test-fast test test-cov clean catalog-init catalog-build seed-atlas-gse-map fetch-geo-sample-metadata catalog-refresh-release summarize-ewas-db-failures retry-ewas-db-failures agent-context references download-cpgcorpus download-cpgcorpus-gse download-ewas-atlas download-ewas-datahub download-ewas-study download-ewas-family download-manifests download-gencode download-cpg-islands setup-methylgpt download-methylgpt export-cpgpt-static export-ewas-sample-info 7b-status 7b-convert-bg

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
	  'seed-atlas-gse-map  Refresh GSE↔Atlas map from NCBI GEO PubMed IDs' \
	  'fetch-geo-sample-metadata  Fetch pilot GEO SOFT → geo_sample_metadata.parquet' \
	  'catalog-refresh-release  Seed Atlas map + deepmat-data-v1 + phenotype census' \
	  '7b-status      Refresh + print Milestone 7B Hub matrix convert progress' \
	  '7b-convert-bg  Background 7B convert watcher (progress docs + finalize)' \
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

seed-atlas-gse-map:
	@if [[ "$${MBS_SKIP_ATLAS_SEED:-0}" != "1" ]]; then \
	  uv run python scripts/seed_atlas_gse_es_map.py; \
	else \
	  printf 'skip seed-atlas-gse-map (MBS_SKIP_ATLAS_SEED=1)\n'; \
	fi

fetch-geo-sample-metadata:
	uv run python scripts/fetch_geo_sample_metadata.py \
	  --studies-file configs/data/geo_backfill_pilot_gse.txt

catalog-refresh-release: seed-atlas-gse-map
	uv run mbs catalog refresh-release
	uv run mbs catalog validate-release
	uv run mbs catalog phenotype-census
	uv run mbs catalog trait-eligibility

summarize-ewas-db-failures:
	uv run python scripts/summarize_ewas_db_download_failures.py

retry-ewas-db-failures:
	bash scripts/retry_ewas_db_download_failures.sh

7b-status:
	bash scripts/status_7b_hub_matrices.sh

7b-convert-bg:
	bash scripts/convert_hub_full_packs_background.sh

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
	uv run python -c "from pathlib import Path; from mbs.paths import DataPaths; from mbs.registry.sample_info import export_family_from_data_root; p=DataPaths.from_environment(); print(export_family_from_data_root(p.data_root, '$(FAMILY)', project_root=p.project_root))"

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
