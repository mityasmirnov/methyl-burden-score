# Static foundation-model features

## Stage 0 principle

Foundation models run offline as feature exporters. The Stage 0 training environment consumes immutable, locus-indexed arrays and does not import CpGPT or MethylGPT at runtime.

This separates:

- expensive external dependencies;
- model and checkpoint licensing;
- GPU feature generation;
- reproducible MBS training.

## CpGPT install and weight download (ops)

CpGPT is an **optional** project extra (not part of default `uv sync --all-groups` / CI). It installs the vendored pin under `vendor/cpgpt` editable.

```bash
cd /data/projects/methyl-burden-score
source scripts/activate_data_environment.sh
uv sync --all-groups --extra cpgpt
```

CpGPT declares `torch<=2.6`. The project keeps a newer MBS torch line via
`[tool.uv] override-dependencies = ["torch>=2.3"]` in `pyproject.toml`.
Linux/Windows torch wheels come from the PyTorch **cu128** index (CUDA 12.8
driver match); see `docs/WORKSPACE.md`.

Download **small** / **human** resources into the project Hugging Face cache
(`$HF_HOME` → `$MBS_CACHE_ROOT/huggingface`; never `$HOME/.cache`):

```bash
uv run --extra cpgpt python - <<'PY'
import os
from cpgpt import download_cpgpt

resources = download_cpgpt(
    model="small",
    species="human",
    cache_dir=os.environ["HF_HOME"],
)
print(resources.checkpoint_path)
print(resources.config_path)
print(resources.dependencies_path)
PY
```

Repeated calls reuse the HF hub cache (~30MB checkpoint + ~5.8G human
dependencies once materialized). Optional: set `HF_TOKEN` for higher Hub rate
limits. Do not commit weights or dependency blobs.

## Default: CpGPT2M sequence-adapter feature

CpGPT begins with a DNA-language-model vector for sequence flanking each CpG and projects it through a methylation-pretrained sequence adapter.

Recommended feature:

```text
source DNA model: NTv2 500M
sequence context: 2001 bp centered on target cytosine
raw dimension: 1024
CpGPT model: small / CpGPT2M
exported dimension: 128
```

Export only the output of `net.encode_sequence()`.

Canonical export command (requires `--extra cpgpt` and the locus registry from
Milestone 2):

```bash
cd /data/projects/methyl-burden-score
source scripts/activate_data_environment.sh
uv sync --all-groups --extra cpgpt
uv run --extra cpgpt mbs features export-cpgpt --feature-set-id cpgpt2m_adapter_128_v1
# or: make export-cpgpt-static
```

Artifacts land under
`$MBS_DATA_ROOT/canonical/static_features/cpgpt2m_adapter_128_v1/`
(`embeddings.zarr`, `loci.parquet`, `artifact.json`). Inspection summary:
`reports/inspection/static_features_cpgpt2m_v1/`.

Do not include:

- sample beta values;
- methylation encoder output;
- CLS/sample embedding;
- transformer token states;
- rotary position encoding;
- task-specific heads.

The adapter output is static, sequence-derived, and available for arbitrary sequence-addressable CpG loci.

**Coordinate convention:** the Methylation Burden Score (MBS) locus registry
stores **1-based** cytosine positions (`GRCh38:chr1:10848`). CpGPT human
dependency keys are **0-based** Ensembl locations (`1:10847`). The exporter
converts automatically.

### CpGPT vs MBS torch pin (why we do not call `CpGPTInferencer`)

**MBS** = Methylation Burden Score (this repo / the `mbs` Python package and CLI).

The Stage 0 training environment pins a newer PyTorch (Linux **cu128**, currently
`torch==2.11.x`) so the GPUs work on the host driver. CpGPT’s optional extra is
installed into the same `.venv`, but:

1. Vendored CpGPT declares `torch<=2.6`; MBS overrides that via
   `override-dependencies` in `pyproject.toml`.
2. Importing `CpGPTInferencer` / `cpgpt.model` pulls
   `torchtune.modules.RotaryPositionalEmbeddings`, which imports
   `torchao.dtypes.nf4tensor` (`NF4Tensor`, `linear_nf4`, …).
3. The torchao wheel resolved with the MBS torch line (**0.18.x** here) no longer
   exposes that NF4 module API, so full CpGPT model import fails with
   `ModuleNotFoundError: torchao.dtypes.nf4tensor`.

Stage 0 only needs the **sequence adapter** (`dna_encoder` /
`encode_sequence()` → 128-d). Export therefore:

- uses `cpgpt.downloads` / Hugging Face cache for checkpoint + human NTv2 deps;
- loads `dna_encoder.*` weights from `small.ckpt` into
  `mbs.static_features.cpgpt_adapter.SequenceAdapterMLP` (local MLPBlock
  equivalent);
- never imports the full CpGPT Lightning module or torchtune at export time.

Training must still **not** import CpGPT; it only reads the frozen
`embeddings.zarr` + `artifact.json`. If upstream CpGPT/torchtune/torchao become
compatible with the MBS torch pin, we can optionally switch back to
`model.net.encode_sequence()` — the artifact contract stays the same.

## Ablation: raw DNA-language-model feature

Export the raw NTv2 vector as a 1024-dimensional baseline. This tests whether CpGPT's methylation pretraining improves the locus representation beyond DNA sequence alone.

## MethylGPT install and weight download (ops)

MethylGPT **cannot** share the main `.venv`. It needs `torchtext`, which only
ships ABI-matched wheels through torch ~2.1–2.4, while the MBS / CpGPT line uses
a newer torch. Use a dedicated project env and keep checkpoints under
`$MBS_DATA_ROOT/raw/methylgpt` (never under `vendor/` or `$HOME`).

```bash
cd /data/projects/methyl-burden-score
source scripts/activate_data_environment.sh
make setup-methylgpt          # creates .venv-methylgpt; pins torch 2.1 + torchtext
make download-methylgpt       # medium checkpoint + type3 probe IDs
```

Equivalent scripts: `scripts/setup_methylgpt_env.sh`,
`scripts/download_methylgpt_weights.sh` (`--medium` default; also `--base`,
`--large`, `--all`).

Layout after download:

```text
$MBS_DATA_ROOT/raw/methylgpt/
├── SOURCE.txt
├── SHA256SUMS                  # local provenance; not required in git
├── pretrained_models/
│   └── methylgpt-medium/       # args.json, *.pt, vocab.json
└── vocab/
    └── probe_ids_type3.csv
```

Import / load smoke check (use the dedicated interpreter):

```bash
.venv-methylgpt/bin/python -c "from methylgpt import MethylGPTModel; print('ok')"
```

The Google Drive medium folder currently ships
`small-best_model_epoch6.pt`; `args.json` confirms `layer_size=128` (medium).
Do not commit weights or `.venv-methylgpt`. Optional: `flash-attn` for
`fast_transformer=True` inference; token-prior export only needs the embedding
table and works with `fast_transformer=False`.

## Ablation: MethylGPT token prior

MethylGPT learns an embedding table indexed by probe vocabulary. This is a population-derived locus identity and co-methylation prior, not a sequence embedding.

Recommended name:

```text
methylgpt_token_prior_medium_128_v1
```

The exporter reads the encoder embedding table, removes special-token rows, maps probe IDs to the canonical GRCh38 locus registry, and records missing or ambiguous mappings.
Load checkpoints from `$MBS_DATA_ROOT/raw/methylgpt/pretrained_models/methylgpt-medium/`
via `.venv-methylgpt` (see ops section above).

Limitations:

- vocabulary-restricted coverage;
- probe-ID rather than coordinate identity;
- possible phenotype enrichment from vocabulary selection;
- no native representation for arbitrary long-read CpGs.

It is therefore an ablation rather than the Stage 0 default.

## Feature artifact layout

```text
static_features/{feature_set_id}/
├── embeddings.zarr
├── loci.parquet
├── normalization.npz          # optional
└── artifact.json
```

`loci.parquet` contains:

```text
embedding_row
locus_id
canonical_key
source_location_key
source_embedding_row
mapping_status
```

## Required manifest fields

```json
{
  "feature_set_id": "cpgpt2m_adapter_128_v1",
  "source_repository": "lucascamillomd/CpGPT",
  "source_commit": "<sha>",
  "model_name": "small",
  "checkpoint_sha256": "<sha256>",
  "configuration_sha256": "<sha256>",
  "source_model": "nucleotide-transformer-v2-500m-multi-species",
  "context_length": 2001,
  "genome_build": "GRCh38",
  "input_dimension": 1024,
  "output_dimension": 128,
  "storage_dtype": "float16",
  "locus_table_sha256": "<sha256>",
  "export_command": "<exact command>",
  "created_at": "<ISO-8601>"
}
```

MethylGPT artifacts additionally require vocabulary and special-token-order hashes.

## Validation

Every export must verify:

- expected row and column counts;
- no NaN or infinite values;
- exact locus ordering;
- checkpoint and configuration compatibility;
- deterministic output under repeated export;
- sampled equality with direct model calls;
- non-degenerate norm and variance distributions;
- explicit missing-locus list;
- genome-build agreement.

## Compression controls

Required feature ablations:

```text
none
raw_ntv2_1024
cpgpt2m_adapter_128
cpgpt2m_adapter_pca32
methylgpt_token_prior
cpgpt_plus_methylgpt
```

PCA is fitted only on static locus vectors and is independent of samples and phenotypes.

## Leakage controls

A static-only model, with methylation values removed, must be evaluated. Strong age or tissue prediction by static features alone indicates platform, coverage, vocabulary-selection, or study leakage.

Compare MethylGPT features against:

- fixed random vectors of equal dimension per locus;
- selected vocabulary with no vector;
- a coverage-matched random locus panel.

## Runtime lookup

Static vectors are stored once per locus and looked up after batch collation:

```python
static = feature_store[batch.locus_row]
cpg_input = torch.cat([batch.sample_features, static, annotations], dim=-1)
```

Do not duplicate static vectors into each sample matrix.