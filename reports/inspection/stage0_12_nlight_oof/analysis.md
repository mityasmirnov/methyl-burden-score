# Milestone 12 — N-light@64 5×6 OOF

Updated: `2026-09-08T17:05+02:00`

- Split: `hub-nine-pack-5fold-v1` (34 234 / 470 studies)
- Product readout: `mbs_enet_nested` (CPU after each GPU job)
- Encoder diagnostic: `mbs_e2e`
- Encoder aux: 9 disease classes + 9 cancer types with nine-pack disease-tissue **n>200**
- Device: GPU 2 (~98 GB), `batch_size=1024` (~87 GB used)
- Budget: **30 epochs**, early-stopping patience **15**

## Progress

| Item | Status |
|------|--------|
| 16-ep plumbing | **archived** `stage0-12-nlight-oof-ep16-f0-r0` |
| 30-ep 5×6 | **running** — `f0-r0` **epoch 19–20 / 30** (1/30 jobs) |
| Aux heads (val) | disease AUROC ~**0.76**; cancer AUROC ~**0.85** (ep 19) |

Runner: `scripts/run_12_nlight_oof.sh` · log `scratch/logs/12_nlight_oof_gpu2_ep30.log`.

Nine-pack 3-fold nested reference (not this 5-fold OOF): **0.368 / 9.88 / 0.803**.

## Live `f0-r0` (30-ep recipe, in flight)

Seed 42. Disease/cancer labels confirmed loaded (9+9). Nested enet not yet
run for this restart.

| Epoch | train_loss | val_loss | val disease AUROC | val cancer AUROC |
|------:|-----------:|---------:|------------------:|-----------------:|
| 1 | 13.54 | 23.20 | 0.41 | 0.57 |
| 10 | 4.88 | 21.11 | 0.58 | 0.83 |
| 15 | 2.44 | 16.74 | 0.73 | 0.85 |
| 16 | 2.36 | **16.59** | 0.75 | 0.84 |
| 19 | 2.01 | 23.30 | **0.76** | **0.85** |

Best val_loss so far: epoch **16**. Training continues toward 30 (patience 15).

## 16-epoch f0-r0 (plumbing, not product)

Seed 42, `max_epochs=16`, patience 5, ATS heads only. Ran 14 epochs,
`best_epoch=9`.

| Readout | Tissue macro-F1 | Age MAE | Sex AUROC |
|---------|----------------:|--------:|----------:|
| `mbs_enet_nested` | **0.300** | **8.63** | **0.880** |
| `mbs_e2e` | 0.286 | 17.74 | 0.748 |

Kept as a 5-fold split / nested-enet wiring check. Remaining 29 jobs would
have inherited 16 epochs from the in-memory runner, so the process was
stopped and retrained under the 30-ep + n>200 recipe.

## Encoder aux heads (locked)

Positives = Hub sidecar `sample_type=disease tissue`. Negatives = pack
`control`. Adjacent-normal and unlabeled nine-pack rows stay **unknown**.
Not `cancer_mask` / `disease_mask`.

**Disease (n>200):** Alzheimer’s 945, schizophrenia 536, SLE 341,
Parkinson’s 333, RA 250, MS 228, psoriasis 211, stroke 204, childhood
asthma 202.

**Cancer types (n>200):** glioma 481, HNSCC 437, AML 437, prostate 401,
stomach 390, breast 343, melanoma 327, HCC 241, sarcoma 239.

Exact-200 types (glioblastoma, ccRCC, cervical, CRC, LUAD, endometrial,
TCC, KIRP, LUSC, thyroid; intellectual disability 200) are out.
