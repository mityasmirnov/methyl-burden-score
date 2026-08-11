# Milestone 6 — hierarchical region model

- **run_id:** `stage0-hier-smoke-maxloci`
- **config:** `configs/experiment/stage0_hier_deeprvat_full.yaml`
- **matrix:** `matrix-hub-age-tissue-sex-full-v1`
- **max_loci:** `2000` (null = full matrix columns)
- **plan:** `docs/plans/milestone-6-hierarchical-region-model.md`

## Topology

- genes: `51`
- regions: `791`
- typed edges: `1634`
- unassigned singleton regions: `434`
- region types: `['promoter_core', 'promoter_proximal', 'five_prime', 'three_prime', 'gene_body', 'unassigned']`
- policy: singleton region_type=unassigned → synthetic gene __unassigned__; Illumina-coordinate-unmapped probes remain matrix-excluded

## Split

- reused flat 5d split: `True`
- train/val/test n: 9489/2074/1985

## External / holdout metrics

- hierarchical: accuracy=0.013256006628003313, mae=22.50695684984937, sex_accuracy=0.49389272437599574
- flat 5d: accuracy=0.6661143330571665, mae=21.976888244485455, sex_accuracy=0.9314922995220393

## Ablations (holdout subset)

- **full:** accuracy=0.0, mae=24.758179014351725, sex_accuracy=0.4360587002096436, n=512
- **drop_unassigned:** accuracy=0.0, mae=24.86836194296756, sex_accuracy=0.4360587002096436, n=512
- **promoters_only:** accuracy=0.0, mae=24.96955034679843, sex_accuracy=0.5828092243186582, n=512
- **gene_body_only:** accuracy=0.0, mae=24.386780095105824, sex_accuracy=0.4171907756813417, n=512

## Artifacts

- metrics: `/data/projects/methyl-burden-score/artifacts/runs/stage0-hier-smoke-maxloci/metrics.json`
- checkpoints: `/data/projects/methyl-burden-score/artifacts/checkpoints/stage0-hier-smoke-maxloci`
- TensorBoard: `/data/projects/methyl-burden-score/artifacts/runs/stage0-hier-smoke-maxloci/tb`
- full uncapped train log: `scratch/logs/hier_full.log`
