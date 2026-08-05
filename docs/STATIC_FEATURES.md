# Static foundation-model features

## Stage 0 principle

Foundation models run offline as feature exporters. The Stage 0 training environment consumes immutable, locus-indexed arrays and does not import CpGPT or MethylGPT at runtime.

This separates:

- expensive external dependencies;
- model and checkpoint licensing;
- GPU feature generation;
- reproducible MBS training.

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

Do not include:

- sample beta values;
- methylation encoder output;
- CLS/sample embedding;
- transformer token states;
- rotary position encoding;
- task-specific heads.

The adapter output is static, sequence-derived, and available for arbitrary sequence-addressable CpG loci.

## Ablation: raw DNA-language-model feature

Export the raw NTv2 vector as a 1024-dimensional baseline. This tests whether CpGPT's methylation pretraining improves the locus representation beyond DNA sequence alone.

## Ablation: MethylGPT token prior

MethylGPT learns an embedding table indexed by probe vocabulary. This is a population-derived locus identity and co-methylation prior, not a sequence embedding.

Recommended name:

```text
methylgpt_token_prior_medium_128_v1
```

The exporter reads the encoder embedding table, removes special-token rows, maps probe IDs to the canonical GRCh38 locus registry, and records missing or ambiguous mappings.

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