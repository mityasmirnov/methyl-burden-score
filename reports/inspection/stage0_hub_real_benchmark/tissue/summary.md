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
