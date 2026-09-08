# Milestone 12 — N-light@64 5×6 OOF

Updated: `2026-09-08T16:20+02:00`

- Split: `hub-nine-pack-5fold-v1` (34 234 / 470 studies)
- Product readout: `mbs_enet_nested` (CPU after each GPU job)
- Encoder diagnostic: `mbs_e2e`
- Encoder aux: 9 disease classes + 9 cancer types with nine-pack disease-tissue **n>200**
- Device: GPU 2 (~98 GB), `batch_size=1024`
- Budget: **30 epochs**, early-stopping patience **15**

## Progress

| Item | Status |
|------|--------|
| 16-ep plumbing | **archived** `stage0-12-nlight-oof-ep16-f0-r0` |
| 30-ep 5×6 | **restarting** (`stage0-12-nlight-oof-f{fold}-r{restart}`) |

Runner: `scripts/run_12_nlight_oof.sh`.

Nine-pack 3-fold nested reference (not this 5-fold OOF): **0.368 / 9.88 / 0.803**.

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
