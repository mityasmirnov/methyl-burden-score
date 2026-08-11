# Stage 0 Milestone 5d — max-N flat DeepRVAT baseline

Public model name: **deepMAT**. Shared flat encoder with **decoupled**
age / tissue / sex phenotype modules and **masked** per-trait loss
(DeepRVAT pattern). Disease/cancer not in scope for 5d.

## Matrices

- `matrix-hub-age-full-v1`: shape `[8374, 482379]` (8374 samples x 482379 loci)
- `matrix-hub-tissue-full-v1`: shape `[5323, 482379]` (5323 samples x 482379 loci)
- `matrix-hub-sex-full-v1`: shape `[2978, 482379]` (2978 samples x 482379 loci)
- **Merged** `matrix-hub-age-tissue-sex-full-v1`: shape `[13548, 482379]` (GSM-union; notes: Multitask merge of matrix-hub-age-full-v1 + matrix-hub-tissue-full-v1 + matrix-hub-sex-full-v1; GSM dedupe=3127)

## Phenotype table

- `sample_phenotype_table_age_tissue_sex_full_v1.parquet`: **13548** samples
- masks: age=10002, tissue=7866, sex=12445

## Train run

- run_id: `stage0-flat-deeprvat-age-tissue-sex-full-v1`
- config: `configs/experiment/stage0_flat_deeprvat_full.yaml`
- split: train=9489 / val=2074 / test=1985 (studies 225 / 60 / 42)
- tissue classes: **47**; genes: **19554**
- best_epoch: `10` (best_val_loss=3.785450618019419)
- external_test tissue accuracy: **0.6661143330571665** (n=1207.0)
- external_test age MAE: **21.976888244485455** years (n=1375.0)
- external_test sex accuracy: **0.9314922995220393** (n=1883.0)

## Artifacts

- run: `$MBS_ARTIFACT_ROOT/runs/stage0-flat-deeprvat-age-tissue-sex-full-v1/`
- checkpoints: `$MBS_ARTIFACT_ROOT/checkpoints/stage0-flat-deeprvat-age-tissue-sex-full-v1/`
- plan: [`docs/plans/milestone-5d-max-n-flat-baseline.md`](../../plans/milestone-5d-max-n-flat-baseline.md)

