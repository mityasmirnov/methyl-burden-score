# Milestone 6 — hierarchical region model

- **run_id:** `stage0-hier-deeprvat-age-tissue-sex-full-v1`
- **config:** `configs/experiment/stage0_hier_deeprvat_full.yaml`
- **matrix:** `matrix-hub-age-tissue-sex-full-v1`
- **max_loci:** `None` (null = full matrix columns)
- **plan:** `docs/plans/milestone-6-hierarchical-region-model.md`

## Topology

- genes: `19554`
- panel (genes + residual slot): `19555`
- regions: `100453`
- typed edges: `440903`
- residual columns: `108070`
- annotation status (study columns): mapped `318965`, multi_mapped `55344`, unmapped `108070`, ambiguous `0`
- region types: `['promoter_core', 'promoter_proximal', 'five_prime', 'three_prime', 'gene_body']`
- policy: mapped loci → typed CpG→region→gene; unmapped/ambiguous → residual path (no `__unassigned__` gene pooling); Illumina-coordinate-unmapped probes retained as residual matrix columns on reconvert

## Split

- reused flat 5d split: `True`
- train/val/test n: 9489/2074/1985

## External / holdout metrics

- hierarchical: accuracy=0.5981772990886496, mae=27.774534330316442, sex_accuracy=0.9341476367498672
- flat 5d: accuracy=0.6661143330571665, mae=21.976888244485455, sex_accuracy=0.9314922995220393
- vs_flat delta: `{'mae': 5.797646085830987, 'accuracy': -0.06793703396851691, 'sex_accuracy': 0.0026553372278279586, 'loss': -0.7721252325626597}`

## Annotation slices (mapped vs residual)

Holdout subset (`ablation_max_samples=512`):

- **full:** accuracy=0.9022, mae=33.89, sex_accuracy=0.9853, n=512
- **mapped_only:** accuracy=0.8978, mae=32.77, sex_accuracy=0.9853, n=512
- **residual_only:** accuracy=0.0, mae=27.45, sex_accuracy=0.4172, n=512

## Ablations (holdout subset)

- **full:** accuracy=0.9022, mae=33.89, sex_accuracy=0.9853, n=512
- **mapped_only:** accuracy=0.8978, mae=32.77, sex_accuracy=0.9853, n=512
- **residual_only:** accuracy=0.0, mae=27.45, sex_accuracy=0.4172, n=512
- **promoters_only:** accuracy=0.8578, mae=51.19, sex_accuracy=0.9790, n=512
- **gene_body_only:** accuracy=0.0711, mae=31.77, sex_accuracy=0.9832, n=512

## Analysis

Uncapped hierarchical train completed (best epoch **13**, early-stopped after 21
epochs; best val loss 3.30). Same 5d folds as flat deepMAT
(`reused_flat_split=true`).

**Vs flat (full external test, n=1985):** hierarchical underperforms flat on
tissue accuracy (−6.8 pp: 0.598 vs 0.666) and age MAE (+5.8 y: 27.8 vs 22.0).
Sex accuracy is essentially tied (+0.3 pp: 0.934 vs 0.931). Multitask loss is
slightly better for hierarchical (−0.77). So Stage 0 reference for phenotype
heads remains flat 5d; hierarchical is the first typed-region + residual
baseline, not a win on holdout accuracy/MAE.

**Mapped vs residual:** on the ablation subset, `mapped_only` ≈ `full` for
tissue/sex, while `residual_only` collapses (tissue acc 0, sex ~chance). The
gene/region path carries almost all predictive signal; the residual path is
retained for completeness and future work, not as a standalone phenotype
encoder in this run.

**Role ablations:** `promoters_only` keeps most tissue signal (0.858) but
hurts age MAE badly (51 y). `gene_body_only` nearly removes tissue accuracy
(0.071) while sex stays high — consistent with sex signal being broadly
distributed / easy, and tissue relying more on promoter-proximal structure in
this setup.

**Caveat:** existing Hub matrix predates Illumina-unmapped residual columns;
the 108 070 residual cols here are gene-unassigned / orphan loci on the current
matrix. Reconvert will add Illumina-coordinate residual probes under the same
path.

## Artifacts

- metrics: `/data/projects/methyl-burden-score/artifacts/runs/stage0-hier-deeprvat-age-tissue-sex-full-v1/metrics.json`
- checkpoints: `/data/projects/methyl-burden-score/artifacts/checkpoints/stage0-hier-deeprvat-age-tissue-sex-full-v1`
- TensorBoard: `/data/projects/methyl-burden-score/artifacts/runs/stage0-hier-deeprvat-age-tissue-sex-full-v1/tb`
- full uncapped train log: `scratch/logs/hier_full.log`
