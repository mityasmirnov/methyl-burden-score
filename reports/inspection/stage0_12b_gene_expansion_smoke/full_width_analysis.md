# Full gene-linked width + CpGPT — first ~20k-gene result

Updated: `2026-09-28` — single run, fold 0. `stage0-12b-cpgpt-full-width-smoke-f0`.

First run at the **~20k-gene product scale** the programme has been aiming at:
`max_loci: null` → **374 309 gene-linked CpG columns / 440 903 edges /
19 554 genes** (vs 2 646 genes at the 65k prefix). CpGPT2M on. Regularization
raised vs the 65k configs (dropout 0.3, weight_decay 1e-3) to counter the
overfit seen in the earlier no-CpGPT full-width attempt.

Did **not** OOM: the VRAM probe calibrated batch **222 → 55**.

## Result (nested enet = product readout)

| arm | genes | nested tissue F1 | nested age MAE | nested sex AUROC |
|---|---:|---:|---:|---:|
| 65k, no CpGPT (n=1) | 2 646 | 0.376 | 9.567 | 0.814 |
| 65k, + CpGPT (n=6 mean) | 2 646 | 0.355 | **8.096** | 0.879 |
| **full width, + CpGPT (n=1)** | **19 554** | **0.388** | 11.671 | **0.962** |

Δ full-width vs 65k+CpGPT: tissue **+0.033**, age **+3.58 (worse)**, sex **+0.083**.

Also on this run: `mbs_e2e` tissue 0.335 / age 16.298 / sex 0.977;
`mbs_linear_probe` tissue 0.359 / age 8.215 / sex 0.977.

## Reading it

**Tissue and sex improve substantially at full gene scale.** Nested tissue
**0.388** matches the long-standing tissue leader in this programme — classical
`C-mvalue-enet-G` at **0.388 ± 0.018** (`TODO_PIPELINE` §"Trustworthy ATS
numbers"), which no neural MBS arm had matched. Sex (0.962) is the highest
recorded on any arm here. This is the first direct support for the
**full-scale hypothesis**: more genes (and eventually non-coding CpGs) buys
representation quality that the 65k prefix cannot.

**Age regresses badly** (8.10 → 11.67) and that is the headline caveat, because
age is the designated primary metric. Two plausible causes, not separated here:

1. **Undertrained / overfit.** val_loss rose 13.9 → 36.4 and early stopping
   fired at epoch 9/16 with `best_epoch` **4**. The model barely trained. The
   dropout/weight-decay bump softened the earlier collapse (71.97 → 36.4) but
   did not fix it.
2. **Nested enet at 19 554 features.** `milestone-12b-full-gene-panel.md` §3
   explicitly warns that nested enet on ~20k MBS columns overfits without
   adequate inner-CV regularization, citing P2-G's nested RBS age blowing up at
   ~13k region features. This run is squarely in that regime, and age is the
   metric that failure mode hits first.

Note the pattern is *internally consistent* with (2): a classification metric
with many weakly-informative features (tissue, sex) benefits from the extra
columns, while a regression metric degrades — which is what an
under-regularized high-dimensional linear readout does.

## Status: promising, not yet a result to build on

- **n=1**, single seed. Measured single-seed spread at the 65k prefix is ±0.19
  age MAE / ±0.010 tissue; at 19 554 features the spread is **unmeasured** and
  plausibly larger.
- Early-stopped at epoch 9 of 16 with best epoch 4 — the comparison is against
  65k arms that trained to epochs 9–16, so this is **not** a matched-budget
  comparison.
- Fold 0 only; N-light only (gene-linked, so **still no non-coding CpGs** — that
  needs cascade's orphan/direct path).

## What would make this conclusive

1. Fix the overfit first: longer patience, stronger regularization, or fewer
   epochs at higher LR — the run needs to actually converge before its metrics
   mean anything.
2. Separate the two causes of the age regression by scoring the **same**
   checkpoint's MBS through a *more strongly regularized* nested fit; if age
   recovers, cause (2) dominates and the encoder is fine.
3. Multi-restart once a converging recipe exists.
