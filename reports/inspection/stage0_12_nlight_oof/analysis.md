# Milestone 12 — N-light@64 5×6 OOF

Updated: `2026-09-10T10:46+02:00`

- Split: `hub-nine-pack-5fold-v1` (34 234 / 470 studies)
- Product readout: `mbs_enet_nested` (CPU after each GPU job; overlapped)
- Encoder diagnostic: `mbs_e2e`
- Encoder aux: 9 disease + 9 cancer types (nine-pack disease-tissue **n>200**)
- Device: GPU 2 (~84–87 GB reserved, `batch_size=1024`)
- Budget: **30 epochs**, early-stopping patience **15**
- **Scope:** HM450 **65k-prefix gene-linked** (**2 646** genes / 51 375 CpGs).
  Not the ~20k-gene product. See
  [`milestone-12b-full-gene-panel.md`](../../../docs/plans/milestone-12b-full-gene-panel.md).

Machine-readable interim: [`summary_interim_24of30.json`](summary_interim_24of30.json).

## Progress

| Item | Status |
|------|--------|
| 16-ep plumbing | **archived** `stage0-12-nlight-oof-ep16-f0-r0` |
| 30-ep 5×6 | **24/30** complete (folds 0–3 × 6); **f4-r0** in flight (~epoch 5/30) |
| Completeness gate | `check_12_oof_completeness.py` → **INCOMPLETE** (6 missing) |
| ETA | ~8 h wall at ~1.34 h/job → ~18:30–19:00 local if uninterrupted |

Runner: `scripts/run_12_nlight_oof.sh` · log `scratch/logs/12_nlight_oof_gpu2_current.log`.

## Interim product readout (24/30) — `mbs_enet_nested`

Grand mean ± population SD over 24 fold×restart jobs (not yet fold-averaged
then mean — same as campaign 3-fold tables for comparability):

| Readout | Tissue macro-F1 | Age MAE | Sex AUROC |
|---------|----------------:|--------:|----------:|
| **`mbs_enet_nested` (24/30)** | **0.325 ± 0.031** | **9.79 ± 0.92** | **0.831 ± 0.057** |
| `mbs_e2e` (diagnostic) | 0.319 ± 0.042 | 15.17 ± 2.79 | 0.814 ± 0.094 |
| Ref: nine-pack **3-fold** nested | 0.368 | 9.88 | 0.803 |

**Takeaways (interim, not final):**

- Nested enet still beats e2e on **age** (~5–6 MAE points) — same pattern as
  the 3-fold reference.
- Tissue F1 is **~0.04 below** the 3-fold nested reference (0.325 vs 0.368).
  Fold structure explains most of it: fold 1 is strong (~0.37); folds 0/3 sit
  ~0.30. Do **not** call this a failure until fold 4 lands and the gate is green.
- Sex holds or improves vs the 3-fold ref.
- Val aux at best epoch (encoder heads): disease AUROC ~**0.80**, cancer
  ~**0.87** (mean over 24 best checkpoints).

### Per-fold restart means (nested)

| Fold | n | Tissue F1 | Age MAE | Sex AUROC |
|-----:|--:|----------:|--------:|----------:|
| 0 | 6 | 0.307 | 9.62 | 0.894 |
| 1 | 6 | 0.372 | 10.08 | 0.763 |
| 2 | 6 | 0.324 | 8.97 | 0.852 |
| 3 | 6 | 0.297 | 10.48 | 0.814 |
| 4 | 0 | — | — | — |

## Encoder aux heads (locked)

Positives = Hub sidecar `sample_type=disease tissue`. Negatives = pack
`control`. Adjacent-normal and unlabeled nine-pack rows stay **unknown**.

**Disease (n>200):** Alzheimer’s 945, schizophrenia 536, SLE 341,
Parkinson’s 333, RA 250, MS 228, psoriasis 211, stroke 204, childhood
asthma 202.

**Cancer types (n>200):** glioma 481, HNSCC 437, AML 437, prostate 401,
stomach 390, breast 343, melanoma 327, HCC 241, sarcoma 239.

## 16-epoch f0-r0 (plumbing, not product)

| Readout | Tissue macro-F1 | Age MAE | Sex AUROC |
|---------|----------------:|--------:|----------:|
| `mbs_enet_nested` | **0.300** | **8.63** | **0.880** |
| `mbs_e2e` | 0.286 | 17.74 | 0.748 |

## Next

1. Finish f4-r0…r5 (~6 jobs).
2. Run `uv run python scripts/check_12_oof_completeness.py --run-prefix stage0-12-nlight-oof --n-folds 5`.
3. Refresh this report + `TODO_PIPELINE` §12 Results/Verdict.
4. **GATE G1–G3** before cascade (G4/10e already FAIL → native P2-G).
