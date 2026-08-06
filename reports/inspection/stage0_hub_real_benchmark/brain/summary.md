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
