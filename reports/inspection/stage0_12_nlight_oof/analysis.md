# Milestone 12 — N-light@64 5×6 OOF

Updated: `2026-09-24T00:00+02:00`

- Split: `hub-nine-pack-5fold-v1` (34 234 / 470 studies)
- Product readout: `mbs_enet_nested` (CPU after each GPU job; overlapped)
- Encoder diagnostic: `mbs_e2e`
- Encoder aux: 9 disease + 9 cancer types (nine-pack disease-tissue **n>200**)
- Device: GPU 2 (~84–87 GB reserved, `batch_size=1024`)
- Budget: **30 epochs**, early-stopping patience **15**
- **Scope:** HM450 **65k-prefix gene-linked** (**2 646** genes / 51 375 CpGs).
  Not the ~20k-gene product. See
  [`milestone-12b-full-gene-panel.md`](../../../docs/plans/milestone-12b-full-gene-panel.md).

Machine-readable final: [`summary.json`](summary.json).

## Progress

| Item | Status |
|------|--------|
| 16-ep plumbing | **archived** `stage0-12-nlight-oof-ep16-f0-r0` |
| 30-ep 5×6 | **30/30 complete** (folds 0–4 × 6 restarts) |
| Completeness gate | `check_12_oof_completeness.py --run-prefix stage0-12-nlight-oof --n-folds 5` → **COMPLETE**, all 30 finite |

Runner: `scripts/run_12_nlight_oof.sh` · log `scratch/logs/12_nlight_oof_gpu2_ep30.log`.

## Final product readout (30/30) — `mbs_enet_nested`

Grand mean ± population SD over 30 fold×restart jobs (not fold-averaged then
meaned — same convention as the 3-fold campaign tables for comparability):

| Readout | Tissue macro-F1 | Age MAE | Sex AUROC |
|---------|----------------:|--------:|----------:|
| **`mbs_enet_nested` (30/30, FINAL)** | **0.316 ± 0.034** | **9.59 ± 0.93** | **0.822 ± 0.054** |
| `mbs_e2e` (diagnostic) | 0.308 ± 0.044 | 15.50 ± 2.86 | 0.791 ± 0.101 |
| Ref: nine-pack **3-fold** nested | 0.368 | 9.88 | 0.803 |

**Takeaways (final):**

- Nested enet beats e2e on **age** (~6 MAE points) and **sex** (~3 pts AUROC)
  — same pattern as the 3-fold reference and the 24/30 interim; confirms
  nested enet as the correct product readout, not `mbs_e2e`.
- Tissue F1 is **0.316**, ~0.05 below the 3-fold nested reference (0.368).
  Fold 4 (the last to land) came in weakest (0.279 nested / 0.263 e2e) and
  pulled the grand mean down from the 24/30 interim figure (0.325). This
  reads as **fold variance**, not a monotonic regression: fold 1 stays
  strong (0.372), fold 2 mid (0.324), folds 0/3/4 sit 0.28–0.31. No sign of
  collapse (no fold near-zero) — call the drop a **soft miss** against the
  reference band, not a failure.
- Age MAE (9.59) and sex AUROC (0.822) both hold at or above the 3-fold ref.
- Val aux at best epoch (encoder heads, from `history[best_epoch]`): cancer
  AUROC **0.858** (mean over all 30 checkpoints); disease AUROC **0.737**
  (mean over **24/30** — fold 3 has no `val_disease_auroc` key in its best
  epoch, i.e. no disease-tissue positives landed in that fold's validation
  slice; not investigated further here, flag for `12c`/data-census follow-up
  if disease aux becomes load-bearing).

### Per-fold restart means (nested)

| Fold | n | Tissue F1 | Age MAE | Sex AUROC |
|-----:|--:|----------:|--------:|----------:|
| 0 | 6 | 0.307 | 9.62 | 0.894 |
| 1 | 6 | 0.372 | 10.08 | 0.763 |
| 2 | 6 | 0.324 | 8.97 | 0.852 |
| 3 | 6 | 0.297 | 10.48 | 0.814 |
| 4 | 6 | 0.279 | 8.76 | 0.789 |

### Per-fold restart means (e2e diagnostic)

| Fold | Tissue F1 | Age MAE | Sex AUROC |
|-----:|----------:|--------:|----------:|
| 0 | 0.320 | 15.38 | 0.882 |
| 1 | 0.380 | 12.61 | 0.758 |
| 2 | 0.289 | 17.48 | 0.754 |
| 3 | 0.286 | 15.23 | 0.864 |
| 4 | 0.263 | 16.78 | 0.698 |

### Per-fold aux AUROC (val, best epoch)

| Fold | Disease AUROC | Cancer AUROC |
|-----:|---------------:|-------------:|
| 0 | 0.774 | 0.825 |
| 1 | 0.885 | 0.897 |
| 2 | 0.746 | 0.934 |
| 3 | — (no positives in val) | 0.815 |
| 4 | 0.542 | 0.820 |

## Encoder aux heads (locked)

Positives = Hub sidecar `sample_type=disease tissue`. Negatives = pack
`control`. Adjacent-normal and unlabeled nine-pack rows stay **unknown**.

**Disease (n>200):** Alzheimer's 945, schizophrenia 536, SLE 341,
Parkinson's 333, RA 250, MS 228, psoriasis 211, stroke 204, childhood
asthma 202.

**Cancer types (n>200):** glioma 481, HNSCC 437, AML 437, prostate 401,
stomach 390, breast 343, melanoma 327, HCC 241, sarcoma 239.

## 16-epoch f0-r0 (plumbing, not product)

| Readout | Tissue macro-F1 | Age MAE | Sex AUROC |
|---------|----------------:|--------:|----------:|
| `mbs_enet_nested` | **0.300** | **8.63** | **0.880** |
| `mbs_e2e` | 0.286 | 17.74 | 0.748 |

## Verdict

**N-light@64 5×6 OOF on the 65k-prefix HM450 panel: DONE.** Nested enet is
the confirmed product readout (age/sex hold or beat the 3-fold reference;
tissue is a soft miss, driven by fold variance rather than collapse). This
closes the N-light arm of Milestone 12 at the **65k-prefix validation**
scope — it does **not** stand in for the ~20k-gene product (12b) or unblock
cascade on its own.

## Next

1. ~~Finish f4-r0…r5~~ **done**.
2. ~~Run `check_12_oof_completeness.py`~~ **done — COMPLETE, 30/30**.
3. ~~Refresh this report + `TODO_PIPELINE` §12 Results/Verdict~~ **done**.
4. **GATE G1–G3** before cascade 5×6 (G4/10e already FAIL → native P2-G is
   the locked topology regardless). G1 = Milestone 11 fold-selected panel /
   trait-universe reselect (`in_progress`); G2 = 12b CpGPT/positional probe
   (`pending`); G3 = 12c platform robustness (`pending`, not started).
