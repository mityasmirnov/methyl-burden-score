# Inspection: GSE35069 EWAS Data Hub → canonical matrix

- matrix_id: `matrix-gse35069-ewasdb-v1`
- platform_id: `HM450`
- processing_level: `gmqn`
- genome_build: `GRCh38`
- shape: `[60, 485470]` (samples x loci)
- dtype: `float32`
- conversion_commit: `a48ab0b8e05db940cc8ea7154691290a29a65aa7`

## Stats

```json
{
  "beta_max": 1.0,
  "beta_mean": 0.5290985703468323,
  "beta_min": 0.0,
  "matrix_paths": {
    "betas": "/data/projects/methyl-burden-score/data/canonical/matrices/matrix-gse35069-ewasdb-v1/betas.zarr",
    "locus_index": "/data/projects/methyl-burden-score/data/canonical/matrices/matrix-gse35069-ewasdb-v1/locus_index.parquet",
    "manifest": "/data/projects/methyl-burden-score/data/canonical/matrices/matrix-gse35069-ewasdb-v1/matrix_manifest.json",
    "root": "/data/projects/methyl-burden-score/data/canonical/matrices/matrix-gse35069-ewasdb-v1",
    "sample_index": "/data/projects/methyl-burden-score/data/canonical/matrices/matrix-gse35069-ewasdb-v1/sample_index.parquet"
  },
  "n_collapsed_probes": 34,
  "n_finite_betas": 28965328,
  "n_mapped_probes": 485504,
  "n_missing_cells": 162872,
  "n_observed_probes": 485512,
  "n_out_of_range": 0,
  "n_samples": 60,
  "n_study_loci": 485470,
  "n_unmapped_probes": 8,
  "roundtrip": {
    "max_abs_diff": 0.0,
    "n_compared": 6144,
    "n_mismatch": 0,
    "ok": true,
    "sample_ids": [
      "GSM861635",
      "GSM861665",
      "GSM861694"
    ]
  },
  "unmapped_probe_ids_head": [
    "cg01759674",
    "cg03455418",
    "cg06102775",
    "cg06677538",
    "cg08417382",
    "cg11598976",
    "cg22940546",
    "cg25005368"
  ]
}
```

## Round-trip

PASS — compared 6144 cells across 3 samples; max_abs_diff=0.0.
