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
