# 7H Track B.4 — ATS cascade pooling seed-43

Generated: `2026-09-08T00:48:18.420993+00:00`

Seed **43** on `hub-ats-7e-3fold-v1` (exercises the seed-offset fix).

| Arm | folds | Tissue F1 | Age MAE | Sex AUROC |
|-----|------:|----------:|--------:|----------:|
| P2-G-s2 | 3 | 0.379 (±0.054) | 19.251 | 0.765 |
| scalar-mean-max-s2 | 3 | 0.352 (±0.055) | 22.200 | 0.706 |
| scalar-max-mean-s2 | 3 | 0.353 (±0.046) | 20.824 | 0.703 |
| vector-mean-max-s2 | 3 | 0.342 (±0.055) | 21.735 | 0.658 |

Compare to seed-42 matched 16-ep screen in `stage0_7g_gene_only_probe/analysis.md`.

