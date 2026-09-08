# Nine-pack P2-G / m-only reference snapshot

Generated: `2026-09-08T09:46:13.872287+00:00`

Auto-exported pair only. Curated campaign board with pooling grid,
ablations, and interpretations: [`analysis.md`](analysis.md).

- Matrix: `matrix-hub-nine-pack-virtual-v1` (virtual multi-store)
- Split: `hub-nine-pack-3fold-v1` (**34234** samples)
- Platform: **HM450 only** — no cross-platform claim
- Phase: `full_reference`

## Loader self-check

`{"matrix_id": "matrix-hub-nine-pack-virtual-v1", "shape": [34234, 482379], "checked_rows": 8, "checked_cols": 64, "finite_fraction": 0.91796875, "packs_in_check": {"matrix-hub-blood-full-v1": 8}, "platform_claim": "HM450_only_no_cross_platform"}`

## Results (`mbs_e2e`, outer test)

| Arm | folds | Tissue F1 | Age MAE | Sex AUROC |
|-----|------:|----------:|--------:|----------:|
| P2-G cascade | 3 | 0.355 (±0.054) | 13.431 (±2.319) | 0.853 (±0.087) |
| one-hop m-only | 3 | 0.273 (±0.061) | 15.057 (±2.634) | 0.742 (±0.135) |

### Fold detail

```json
{
  "P2-G": [
    {
      "eval_split": "test",
      "tissue_f1": 0.36363491876843274,
      "age_mae": 12.796308918815143,
      "sex_auroc": 0.9496811792587474,
      "best_epoch": 12,
      "fold": 0,
      "run_id": "stage0-7h-nine-pack-P2-G"
    },
    {
      "eval_split": "test",
      "tissue_f1": 0.29745126516793496,
      "age_mae": 16.00099787320679,
      "sex_auroc": 0.8297151660378583,
      "best_epoch": 4,
      "fold": 1,
      "run_id": "stage0-7h-nine-pack-P2-G"
    },
    {
      "eval_split": "test",
      "tissue_f1": 0.4039206479276686,
      "age_mae": 11.494403423776612,
      "sex_auroc": 0.7795242257598387,
      "best_epoch": 10,
      "fold": 2,
      "run_id": "stage0-7h-nine-pack-P2-G"
    }
  ],
  "m_only": [
    {
      "eval_split": "test",
      "tissue_f1": 0.342626863679713,
      "age_mae": 12.072589709898208,
      "sex_auroc": 0.8901185613370428,
      "best_epoch": 13,
      "fold": 0,
      "run_id": "stage0-7h-nine-pack-m-only-f0"
    },
    {
      "eval_split": "test",
      "tissue_f1": 0.24524517234989326,
      "age_mae": 17.05477789444737,
      "sex_auroc": 0.7083240108010432,
      "best_epoch": 6,
      "fold": 1,
      "run_id": "stage0-7h-nine-pack-m-only-f1"
    },
    {
      "eval_split": "test",
      "tissue_f1": 0.23045663488369839,
      "age_mae": 16.044592541373653,
      "sex_auroc": 0.6262277663475293,
      "best_epoch": 7,
      "fold": 2,
      "run_id": "stage0-7h-nine-pack-m-only-f2"
    }
  ]
}
```

## Notes

- Full cohort (no subsample).
- 7H remains open after this smoke; trait expansion still gated on ≥1k-per-arm census.
- Stage B / Milestone 7 OOF stay blocked.
