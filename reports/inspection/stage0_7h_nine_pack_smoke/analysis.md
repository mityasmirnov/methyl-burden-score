# 7H Phase 3 — nine-pack smoke / reference

Generated: `2026-09-07T10:15:06.371671+00:00`

- Matrix: `matrix-hub-nine-pack-virtual-v1` (virtual multi-store)
- Split: `hub-nine-pack-3fold-v1` (**34234** samples)
- Platform: **HM450 only** — no cross-platform claim
- Phase: `smoke`

## Loader self-check

`{"matrix_id": "matrix-hub-nine-pack-virtual-v1", "shape": [34234, 482379], "checked_rows": 8, "checked_cols": 64, "finite_fraction": 0.91796875, "packs_in_check": {"matrix-hub-blood-full-v1": 8}, "platform_claim": "HM450_only_no_cross_platform"}`

## Results (`mbs_e2e`, outer test)

| Arm | folds | Tissue F1 | Age MAE | Sex AUROC |
|-----|------:|----------:|--------:|----------:|
| P2-G cascade | 1 | 0.343 | 15.038 | 0.891 |
| one-hop m-only | 1 | 0.119 | 18.060 | 0.735 |

### Fold detail

```json
{
  "P2-G": [
    {
      "eval_split": "test",
      "tissue_f1": 0.3433791888216678,
      "age_mae": 15.037650264341393,
      "sex_auroc": 0.8914664301304904,
      "best_epoch": 3,
      "fold": 0,
      "run_id": "stage0-7h-nine-pack-P2-G-smoke"
    }
  ],
  "m_only": [
    {
      "eval_split": "test",
      "tissue_f1": 0.11920655025403026,
      "age_mae": 18.060063433191566,
      "sex_auroc": 0.7352226095439334,
      "best_epoch": 3,
      "fold": 0
    }
  ]
}
```

## Notes

- Full cohort (no subsample).
- 7H remains open after this smoke; trait expansion still gated on ≥1k-per-arm census.
- Stage B / Milestone 7 OOF stay blocked.
- Smoke budgets are **3 epochs only** — not comparable to ATS 15/16-ep refs.
- See also [`trait_adequacy.md`](trait_adequacy.md) (Track B.5).
- Full 3-fold P2-G + m-only: `scripts/run_7h_nine_pack_smoke.py --phase full`
  (log `scratch/logs/7h_nine_pack_full.log`).
