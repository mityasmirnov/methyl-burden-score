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
