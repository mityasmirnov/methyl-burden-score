# GATE G2 — CpGPT2M static-embedding ablation (single fold, single seed)

Updated: `2026-09-24`

- Question: do precomputed CpGPT2M sequence-adapter embeddings (128-dim,
  `data/canonical/static_features/cpgpt2m_adapter_128_v1`) concatenated onto
  N-light `flat_region`/`m_only` improve tissue/age/sex vs the matched
  no-CpGPT baseline?
- Split `hub-nine-pack-3fold-v1`, fold 0 only, 16 epochs, 65 536-locus prefix
  (2 646 genes / 51 375 CpG columns). Isolated ablation: the only config diff
  between the two runs is `stage0.use_cpgpt_static_features` /
  `features.static_feature_set` — `batch_size: auto` in the CpGPT run
  auto-calibrated to **256**, matching the baseline's fixed `batch_size: 256`,
  so batch size is not a confound despite the nominal config difference.
- Baseline: `stage0-7h-nine-pack-m-only-wide-f0`. CpGPT arm:
  `stage0-12b-cpgpt-nlight-smoke-f0`. `mbs_enet_nested` for the CpGPT arm was
  not computed inline (config had `stage_a_include_mbs_enet: false`); fit
  post-hoc via `scripts/eval_mbs_enet_from_scores.py --nested` on the saved
  scores for a fair readout.

## Result

| Readout | Arm | Tissue macro-F1 | Age MAE | Sex AUROC |
|---|---|---:|---:|---:|
| `mbs_e2e` | baseline (no CpGPT) | 0.335 | 16.66 | 0.888 |
| `mbs_e2e` | + CpGPT2M | 0.341 | **10.40** | **0.967** |
| `mbs_enet_nested` | baseline (no CpGPT) | **0.376** | 9.57 | 0.814 |
| `mbs_enet_nested` | + CpGPT2M | 0.340 | **8.44** | **0.904** |

**Reading it:** the `mbs_e2e` deltas (age −6.3 MAE, sex +0.08 AUROC) are
large enough to look almost too good for a single static-feature addition —
that's exactly why `mbs_e2e` isn't the gate's decision metric. Under
`mbs_enet_nested` (this project's product readout), the picture is more
modest and mixed: **age MAE improves** (9.57→8.44, ~12% better), **sex AUROC
improves** (0.814→0.904), but **tissue macro-F1 regresses slightly**
(0.376→0.340, ~−0.036). Per the G2 gate criterion ("promote only if
nested/e2e wins"), this is a **partial win**, not a clean one: 2 of 3
nested-enet metrics improve, tissue does not.

## Caveats (do not over-read this)

- **Single fold, single seed.** No restart/fold variance estimate exists yet
  — the 30/30 N-light 5×6 campaign showed per-fold nested-tissue swinging
  ±0.05 from restart/fold alone (`stage0_12_nlight_oof/analysis.md`), which
  is the same order of magnitude as this tissue regression. Cannot yet rule
  out this being noise rather than a real CpGPT-driven tissue cost.
- Not yet wired into the P2-G cascade encoder (`cascade_loop.py` has no
  static-feature plumbing) — this result is N-light-only.
- 65k-prefix scope, same as the rest of Milestone 12 — not the ~20k-gene
  product panel.

## Verdict

**Not a gate pass or fail yet** — genuinely promising on age/sex, open
question on tissue. Before promoting CpGPT into the N-light product or
wiring it into cascade, this needs at minimum a multi-restart re-run on
fold 0 (matching the N-light 5×6 restart protocol) to separate a real
tissue cost from restart noise, ideally a second fold. Recommend queuing
that as the next GPU0 job rather than treating this single run as decisive.

## User decision (2026-09-24)

User's read: the age/sex gains are compelling enough to want CpGPT static
embeddings as a **main-architecture direction going forward** (N-light now,
cascade later), not just an optional ablation. Recorded as the stated
direction.

**Age is the primary decision metric for this call** (per user direction;
consistent with existing project vocabulary — age-primary framing already
appears in `docs/plans/milestone-7g-prime-age-seed-mask.md`/Milestone 9c,
and `mbs_enet_nested`'s age win is called out project-wide as "the dramatic
age win" in `docs/plans/milestone-10e-staged-rbs-mbs-training.md`). Under
that weighting, the nested age MAE improvement (9.57→8.44) is the headline
result and the tissue regression (−0.036) is a secondary concern, not a
tied co-equal metric — it does not by itself block treating CpGPT as the
leading direction.

This does **not**, however, waive the multi-restart condition already
stated: single fold/seed is still not sufficient basis to lock in an
architecture change on its own, consistent with how this project already
treated the P2-G vector warm-start (fold-0-promising, then failed the fair
multi-run re-test — see `docs/ARCHITECTURE_BENCHMARKS.md`). The
multi-restart re-run (6 seeds, fold 0, `stage0-12b-cpgpt-nlight-smoke-f0-r0..5`)
is running now on GPU0; see
[`../stage0_12b_cpgpt_multirestart/summary.json`](../stage0_12b_cpgpt_multirestart/summary.json)
once complete. Next after that holds up: cascade `cascade_loop.py`
static-feature plumbing.
