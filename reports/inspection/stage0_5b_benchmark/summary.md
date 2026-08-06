# Stage 0 Milestone 5b benchmark report
Public model name: **deepMAT** (CLI/package remain `mbs` / `methyl-burden-score`).
## Registry and sample-info
- Registry: `configs/data/phenotype_registry.yaml`
- Checksums: `$MBS_DATA_ROOT/canonical/registries/download_checksums.parquet`
- **age**: 8374 samples, 143 studies, platforms={'450K': 8374}
- **tissue**: 5323 samples, 258 studies, platforms={'450K': 5323}
- **disease**: 14501 samples, 209 studies, platforms={'450K': 14501}

## Runs
| Run | Split | Final val acc/loss | TensorBoard |
|-----|-------|--------------------|-------------|
| `stage0-flat-gse35069-v1` | None | acc=0.1 loss=2.315625 | False |
| `stage0-5b-tissue-holdout-fixture` | study_grouped | acc=1.0 loss=0.006747476290911436 | True |
| `stage0-5b-age-holdout-fixture` | study_grouped | acc=1.0 loss=0.006747476290911436 | True |

## Study-holdout fixture studies
- `stage0-5b-tissue-holdout-fixture`: train=['GSE_A'] val=['GSE_B'] external_test=['GSE_C']
- `stage0-5b-age-holdout-fixture`: train=['GSE_A'] val=['GSE_B'] external_test=['GSE_C']

## Notes
- GSE35069 remains the real pilot smoke benchmark (donor-grouped).
- Age/tissue study-holdout fixtures exercise study-grouped splits + TensorBoard; full Hub profile matrix convert is documented as follow-on (packs already on disk).
- Disease sample-info exported; disease profile zip not required for 5b gate when subset notes are recorded.
- EWAS Atlas remains validation-only.

## Metrics contract
- age: MAE, RMSE (`mbs.evaluation.regression_metrics`)
- binary: AUROC, AUPRC
- multiclass: macro-F1, balanced accuracy, confusion matrix
- always report by holdout study / platform for multi-study runs
