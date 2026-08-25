# 7C Hub disease multilabel smoke

- Run id: `stage0-flat-hub-disease-multilabel-smoke-v1`
- Matrix: `matrix-hub-disease-full-v1` (`max_loci=256`, `max_samples=48`, 1 epoch)
- Task: multitask / masked disease BCE
- Disease labels in vocab: 30
- Holdout loss ≈ 0.70 (BCE active; labeled samples preferred under `max_samples`)
- Holdout `disease_auroc` / `disease_auprc`: often absent under ADR unknown≠control
  (observed cells are mostly positives; binary AUROC needs both classes)
- Trainer still emits `auroc`/`auprc`/`ece` for binary sex/tissue when both classes
  appear; `disease_auroc` when any label has both classes under the observation mask
- Score polarity: see `score_manifest.json`

Config: `configs/experiment/stage0_flat_hub_disease_multilabel_smoke.yaml`
