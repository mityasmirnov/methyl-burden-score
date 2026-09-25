# GATE G2 — CpGPT 6-restart confirmation (fold 0)

Updated: `2026-09-25` — **complete, 6/6 restarts** (seeds 42–47).

Purpose: the single-run G2 probe showed an age/sex win with a tissue
regression, but with n=1 it could not separate signal from restart noise. This
run repeats the *same* config 6× on fold 0 to measure restart spread.

- Config `stage0_12b_cpgpt_nlight_smoke.yaml` (CpGPT2M on), fold 0 of
  `hub-nine-pack-3fold-v1`, 65 536-locus prefix, 16 epochs.
- Baseline for comparison: `stage0-7h-nine-pack-m-only-wide-f0` (no CpGPT),
  **single run** — see caveat below.
- Machine-readable: [`summary.json`](summary.json).

## Per-restart results

| restart | best epoch | nested tissue | nested age MAE | nested sex | e2e tissue | e2e age MAE | e2e sex |
|--------:|-----------:|--------------:|---------------:|-----------:|-----------:|------------:|--------:|
| r0 | 5  | 0.361 | 8.011 | 0.852 | 0.248 | 13.636 | 0.867 |
| r1 | 12 | 0.371 | 7.964 | 0.878 | 0.342 |  9.684 | 0.957 |
| r2 | 16 | 0.339 | 8.359 | 0.892 | 0.339 |  9.760 | 0.968 |
| r3 | 15 | 0.349 | 8.140 | 0.873 | 0.334 | 13.592 | 0.964 |
| r4 | 10 | 0.356 | 7.955 | 0.877 | 0.309 | 11.993 | 0.935 |
| r5 | 15 | 0.355 | 8.148 | 0.902 | 0.321 | 10.811 | 0.966 |

| Readout | mean ± sd | min–max |
|---|---|---|
| **nested tissue F1** | **0.355 ± 0.010** | 0.339–0.371 |
| **nested age MAE** | **8.096 ± 0.140** | 7.955–8.359 |
| **nested sex AUROC** | **0.879 ± 0.015** | 0.852–0.902 |
| e2e tissue F1 | 0.316 ± 0.032 | 0.248–0.342 |
| e2e age MAE | 11.579 ± 1.630 | 9.684–13.636 |
| e2e sex AUROC | 0.943 ± 0.036 | 0.867–0.968 |

## Finding 1 — `mbs_e2e` is ~11× noisier than the product readout

Age MAE restart sd: **1.630 (e2e) vs 0.140 (nested)**. Tissue: 0.032 vs 0.010.

This has a direct methodological consequence, and it invalidated an experiment
design of ours before it ran: the CpGPT architecture sweep was originally
ranking 6 single-seed variants on **e2e** age MAE. At sd ≈ 1.6, a variant needs
a ~3 MAE gap to be distinguishable, so that ranking would have been mostly
noise and could have promoted a lucky seed. The sweep now ranks on nested
(commit `9bd1204`). **Do not rank single-seed configs on `mbs_e2e`.**

The two readouts also disagree in *level*, not just spread — e2e flatters sex
(0.943 vs 0.879) and penalises tissue (0.316 vs 0.355) relative to nested. They
are not interchangeable proxies.

## Finding 2 — nested is stable even though checkpoint selection is not

`best_epoch` ranges 5–16 across restarts, yet nested varies by only ±0.14 MAE.
So the nested elastic-net on frozen scores largely absorbs checkpoint-selection
noise. This partly defuses the concern flagged in
`scripts/check_12_oof_completeness.py` (that `primary_evaluation:
mbs_enet_nested` describes the *reported* readout while checkpoint selection
actually uses `validation_tissue_macro_f1_then_age_mae`): for the nested
readout, which epoch got picked matters much less than expected.

## Finding 3 — CpGPT verdict (nested)

| Metric | baseline (n=1) | CpGPT (n=6) | Δ | baseline vs CpGPT range |
|---|---:|---:|---:|---|
| age MAE | 9.567 | **8.096 ± 0.140** | **−1.47** | far outside (worse than all 6) |
| sex AUROC | 0.814 | **0.879 ± 0.015** | **+0.065** | far outside (worse than all 6) |
| tissue F1 | **0.376** | 0.355 ± 0.010 | **−0.021** | outside (better than all 6) |

- **Age: large, reliable win** (~15% lower MAE; the baseline is nowhere near the
  CpGPT distribution).
- **Sex: reliable win** (+0.065 AUROC).
- **Tissue: a real regression, not noise.** The baseline (0.376) is better than
  *every* CpGPT restart (max 0.371) and CpGPT's sd is only 0.010. But the
  magnitude is **smaller** than the single-run probe implied (−0.021, not
  −0.036), and the single-run nested age figure (8.44) was also pessimistic
  versus the 6-restart mean (8.096).

**Under the age-primary weighting** (user direction, 2026-09-24), this is a net
win: a decisive age gain plus a sex gain against a small, quantified tissue
cost.

## Caveat that still stands

**The baseline is a single run.** Everything above characterises CpGPT's own
restart distribution and asks whether the baseline point falls inside it. It
does not put error bars on the *baseline*, so the Δ magnitudes are not
confidence intervals. Nested spread is small (sd ≈ 0.14 MAE, 0.010 F1), so a
single baseline is probably a fair point estimate — but the clean fix is a
matched 6-restart baseline (same config, CpGPT off), which is cheap and is the
recommended next confirmation.

Scope unchanged: fold 0 only, 65k-prefix, N-light. Not cascade (CpGPT cascade
plumbing landed separately in `efb8518`/`83d5909`/`f46f0ce` and is unproven at
scale), and not the ~20k-gene product panel.
