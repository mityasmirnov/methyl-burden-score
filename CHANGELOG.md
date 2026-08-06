# Changelog

All notable changes to this project are documented here.

The project follows semantic versioning once the first public API is released.

## Unreleased

### Added

- Stage 0 Milestone 3: offline CpGPT2M sequence-adapter static features
  (`mbs features export-cpgpt` / `make export-cpgpt-static`) writing
  `canonical/static_features/cpgpt2m_adapter_128_v1/` with schema-valid
  `artifact.json`. Export loads `dna_encoder` weights via
  `mbs.static_features.cpgpt_adapter` because full `CpGPTInferencer` import
  fails on the MBS torch / torchtune / torchao pin (documented in
  `docs/STATIC_FEATURES.md`).
- Stage 0 Milestone 4: EWAS Data Hub pilot matrix convert
  (`mbs matrix convert`) for `GSE35069` → `matrix-gse35069-ewasdb-v1`, plus
  `make download-ewas-study` and round-trip inspection report.

### Changed

- Pin Linux/Windows `torch` (and `torchvision`) to the PyTorch **cu128** index
  so wheels match host CUDA 12.8 / driver 570.x. Resolves
  `torch.cuda.is_available() == False` caused by PyPI's default `cu130` wheel.
  Documented in `docs/WORKSPACE.md` and `docs/TODO_PIPELINE.md`.

### Previously added

- Optional `cpgpt` extra (`uv sync --extra cpgpt`) installing the vendored
  CpGPT pin with a torch override so MBS can keep torch ≥2.3; install/download
  ops documented in `docs/STATIC_FEATURES.md` and `docs/TODO_PIPELINE.md`.
- Stage 0 project scope and architecture documentation.
- `/data`-only workspace policy for server development.
- Python package, CLI, SQL catalog, configuration, and testing scaffolds.
- Flat and hierarchical permutation-invariant model primitives.
- Static feature manifests for CpGPT and MethylGPT export workflows.
- Study-grouped cross-fitting and data-leakage design contracts.

### Not yet implemented

- Production CpGCorpus ingestion.
- Canonical GRCh38 annotation graph build.
- Foundation-model feature export.
- Full cross-fitted age and tissue training.
- Episignature and epivariant workflows.
