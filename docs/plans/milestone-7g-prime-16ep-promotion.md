# Plan: Matched 16-epoch promotion screen (7G′ Stage A)

Status: **complete** (2026-09-04, `87b22c3`; per-arm provenance landed
`e644ed4` on 2026-09-05). `promotion_decision.json`: `next_gate:
retain_pooling_2x2`, `recommendation: "Retain full 2×2 pooling result; no
pooling lock. Proceed to age-primary seed-mask."` Age-primary seed-mask CUDA
is now **unblocked and running** — see
[`milestone-7g-prime-age-seed-mask.md`](milestone-7g-prime-age-seed-mask.md).

Parent: [`milestone-7g-prime-stage-a-deeprvat-screen.md`](milestone-7g-prime-stage-a-deeprvat-screen.md).
Report: `reports/inspection/stage0_7g_gene_only_probe/analysis.md`.

## Goal

Match epoch budgets before declaring scalar cascade better than vector or
one-hop. Keep **P2-G as the current reference, not a pooling lock.** Fixed
`rbs_enet` (α=0.1, l1_ratio=0.5, no scaler) is diagnostic only — age collapse
with stable sex on vector RBS is a sparse-penalty failure, not evidence the
vector representation is weak.

## Locked decisions

| Choice | Decision |
|--------|----------|
| Ceiling | 16 epochs; patience 5; tissue-primary checkpoint |
| Loss | Keep tissue-primary (`λ_t=3`, `λ_a=0.3`, `λ_s=1`) |
| Panel | `explicit_only` 51 375 CpGs |
| Microbatch | one-hop batch 256 / 16M tokens; cascade batch **128** |
| References | Keep P2-G / P4-G (15 ep); do not rerun |
| Vector max→max | Do **not** promote yet |
| Nested enet | Post-hoc CPU; train-fold StandardScaler + inner-val α/l1 |

## Queue progress (final, 2026-09-04)

| Arm | Status | Notes |
|-----|--------|-------|
| `N-light-gene-max` | **done** 3/3 | Tissue F1 0.336; age MAE 21.6 — below P2 |
| `N-light-gene-mean` | **done** 3/3 | Tissue F1 **0.378** e2e; age MAE 17.1 — edges out P2-G (0.373) on tissue F1; `one_hop_mean_near_p2` rule fired |
| `N-cascade-scalar-mean-max` | **done** 3/3 | Tissue F1 0.346; age MAE 20.4 |
| `N-cascade-scalar-max-mean` | **done** 3/3 | Tissue F1 0.369; age MAE 20.8 |
| `N-cascade-vector-mean-max` | **done** 3/3 | Tissue F1 0.360; age MAE 21.9 |
| `P2-G` (reference, 15 ep) | reference | Tissue F1 0.373; age MAE 15.6 — best age MAE of the neural set |

No arm beats classical (`C-mvalue-enet-G` F1 0.388; `C-mvalue-ridge-G` age
MAE 6.5) outright. Among neural arms, `N-light-gene-mean` and `P2-G` are
statistically tied on tissue F1; `nothing_beats_p2_or_classical` did not
fire, so the full 2×2 cascade pooling grid stays retained rather than locked
to one pooling choice.

Driver: `scripts/run_7g_16ep_promotion_resume.sh`. **Nested / fixed elastic-net
is post-hoc only** — do not run `eval_mbs_enet_from_scores.py` inside the GPU
queue; refresh the report from e2e + linear/Ridge probes first, then enet
offline.

## Decision rules

- If one-hop max ≈ P2-G across three folds → preferred smaller DeepRVAT-like architecture.
- If one-hop **mean** ≈ P2-G (as in the finished 16-ep mean run) → document as a
  viable smaller topology candidate; do not skip the remaining cascade jobs.
- If vector mean→max improves age/sex e2e but still loses after gene pooling → typed-RBS aggregation, not scalar MBS.
- If scalar mixed pooling closes the gap to P2 → retain the full 2×2 pooling result.
- If nothing beats P2-G or raw-CpG baselines → stop architecture sweeps; start age-primary seed-gene.

After the queue finishes: `write_7g_gene_only_probe_report.py` →
`apply_7g_16ep_decision.py` → inspect `promotion_decision.json`.

## Next steps

1. Let the resume script finish scalar mean→max (folds 1–2), then max→mean, then vector mean→max.
2. Refresh report + `promotion_decision.json` (include light-mean 16-ep near-P2).
3. Optional parallel CPU nested enet on finished runs.
4. On unlock: seed-mask CUDA with `--reuse-panels` (audited fold-0 panel).

## Non-goals

- Launching seed-mask GPU training before the rules fire
- Stage B CpG-panel GPU / Milestone 7
- Rerunning P2/P4 for exact 16-ep equality
- Promoting vector max→max
- Killing the 16-ep PIDs to free GPU 0 early
