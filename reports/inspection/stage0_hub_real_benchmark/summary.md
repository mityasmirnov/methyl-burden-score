# Stage 0 Hub real-matrix benchmark (study-grouped)

Public model name: **deepMAT**. Package/CLI unchanged (`mbs`).

These runs use **EWAS Data Hub profile-pack → canonical matrix** subsets,
not synthetic fixtures. TensorBoard + `metrics.jsonl` are enabled on every run.

# Hub real benchmark — age

- run_id: `stage0-hub-age-studyholdout-v1`
- matrix: `matrix-hub-age-studyholdout-v1`
- studies: `GSE51032, GSE56105, GSE55763, GSE78874`
- task: `regression`
- split: train=['GSE51032', 'GSE56105'] val=['GSE55763'] test=['GSE78874']
- best_epoch: `9`
- best_val_loss: `1.0472601579103087`
- final_val_accuracy: `0.0`
- final_val_mae: `0.9190104024750846`
- external_test: `{"accuracy": 0.0, "loss": 2.532346995600632, "mae": 24.154735084170746, "mae_note": "years (destandardized)", "n_samples": 35}`


# Hub real benchmark — tissue

- run_id: `stage0-hub-tissue-studyholdout-v1`
- matrix: `matrix-hub-tissue-studyholdout-v1`
- studies: `GSE58885, GSE52401, GSE97628, GSE78874, GSE75248`
- task: `multiclass`
- split: train=['GSE52401', 'GSE58885', 'GSE97628'] val=['GSE78874'] test=['GSE75248']
- best_epoch: `1`
- best_val_loss: `2.110267857142857`
- final_val_accuracy: `0.0`
- final_val_mae: `0.0`
- external_test: `{"accuracy": 0.0, "loss": 2.0955357142857145, "mae": 0.0, "mae_note": null, "n_samples": 35}`


# Hub real benchmark — blood

- run_id: `stage0-hub-blood-studyholdout-v1`
- matrix: `matrix-hub-blood-studyholdout-v1`
- studies: `GSE56105, GSE56046, GSE51032, GSE56581, GSE97628`
- task: `multiclass`
- split: train=['GSE51032', 'GSE56046', 'GSE56105'] val=['GSE56581'] test=['GSE97628']
- best_epoch: `1`
- best_val_loss: `1.9292410714285715`
- final_val_accuracy: `0.0`
- final_val_mae: `0.0`
- external_test: `{"accuracy": 0.0, "loss": 1.9212053571428571, "mae": 0.0, "mae_note": null, "n_samples": 35}`


# Hub real benchmark — brain

- run_id: `stage0-hub-brain-studyholdout-v1`
- matrix: `matrix-hub-brain-studyholdout-v1`
- studies: `GSE64509, GSE59457, GSE80970, GSE98203, GSE66351`
- task: `multiclass`
- split: train=['GSE59457', 'GSE64509', 'GSE80970'] val=['GSE98203'] test=['GSE66351']
- best_epoch: `1`
- best_val_loss: `3.577232142857143`
- final_val_accuracy: `0.0`
- final_val_mae: `0.0`
- external_test: `{"accuracy": 0.0, "loss": 3.5973214285714286, "mae": 0.0, "mae_note": null, "n_samples": 35}`


## Label harmonization / design notes

- **blood** `phenotype_value` is per-sample cell-fraction strings — benchmark uses the `tissue` column instead.
- **tissue / blood study-holdout:** each holdout study is a **single tissue type absent from train**. Zero holdout accuracy is expected for closed-set multiclass CE (unseen class IDs). This run validates leakage-safe splits + logging, not biological tissue prediction. Prefer multi-tissue studies or open-set / binary heads before interpreting biology.
- **age** `final_val_mae` in run JSON is in **train-fold standardized** units; `external_test.mae` is destandardized years.
- **disease** profile zip still downloading; age-pack GSM overlap for candidate studies is control-only (no case labels). Convert after pack completes; map empty→`control`; fix `ulcerative colitis` / `Ulcerative colitis` casing.
- **cancer** profile zip incomplete — matrix convert deferred.
