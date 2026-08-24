# Stage 0 / Milestone 7C architecture (fixture evidence)

Fixture-only report. No Hub-scale retraining. Graph id target: `graph-grch38-gencode38-cgi-tile-v2`.

## Sampler

- Deterministic `(seed=7, epoch=1)` → 4 batches
- Token budget respected (singleton overflow allowed): `True`

## Constraint-aware split

- Donor `D` samples `s1`/`s2` same role: `True` (`train`)
- Constraint tallies: `{"age_quantile": {"external_test": {"unknown": 1}, "train": {"unknown": 3}, "validation": {"unknown": 1}}, "case_control": {"external_test": {"unknown": 1}, "train": {"unknown": 3}, "validation": {"unknown": 1}}, "n_split_donors": 4, "platform": {"external_test": {"EPIC": 1}, "train": {"EPIC": 1, "HM450": 2}, "validation": {"EPIC": 1}}, "task_mask": {"external_test": {"age=0,tissue=0,sex=0": 1}, "train": {"age=0,tissue=0,sex=0": 3}, "validation": {"age=0,tissue=0,sex=0": 1}}, "tissue_class": {"external_test": {"brain": 1}, "train": {"blood": 2, "liver": 1}, "validation": {"brain": 1}}}`

## Graph v2 (CGI RBS + CpG-count tiles)

- Region system counts: `{"gene": 1, "rbs": 1, "tbs": 2}`
- Loci per system: `{"gene": 1, "rbs": 1, "tbs": 3}`
- Unmapped locus 6 excluded from tiles/RBS: `True`
- Default tile target CpGs (production): `50`

## Orientation (ADR 0008)

- Synthetic anti-correlated MBS → polarity `flipped`, family `predictive_mbs`

## Direct branch

- Elastic-net loci retained (min_studies=2): `6`

## Matched encoder YAML defaults resolved

```json
{
  "activation": "gelu",
  "cpg_hidden_dim": 64,
  "dropout": 0.1,
  "layer_norm": true
}
```

## Tests

See `tests/unit/test_stage0_7c.py` (sampler, heads, splits, controls, graph-v2, arms, disease masked BCE, matched encoder, orientation train-path, long-form multi-hot).

## Residual polish (orientation + long-form)

- Overfit flat train writes `score_manifest.json` with ADR 0008 polarity from train-fold gene-mean MBS vs signed gene-mean M (`hyper_aligned` or `flipped`).
- `load_longform_multilabel`: multi-hot from repeated `sample_id`; unknown sample → all-False mask (not control).
- Config: `configs/experiment/stage0_flat_hub_disease_multilabel.yaml` (Hub smoke when `matrix-hub-disease-full-v1` exists).

## Still waiting on 7B

- Hub disease/cancer full matrices (`sample_index.parquet` + long-form sidecar) for real multilabel BCE smoke.
- Do not use `matrix-hub-disease-from-agepack-v1` as the multi-label path.
- Other 7C leftovers (AUROC trainer emission, full-genome graph-v2, multi-system hier, true RBS/TBS arms): see `docs/plans/milestone-7c-supervised-architecture.md`.
