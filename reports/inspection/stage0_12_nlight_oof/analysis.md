# Milestone 12 — N-light@64 5×6 OOF

Updated: `2026-09-08T15:50+02:00`

- Split: `hub-nine-pack-5fold-v1` (34 234 / 470 studies)
- Product readout: `mbs_enet_nested` (CPU after each GPU job)
- Encoder diagnostic: `mbs_e2e`
- Device: GPU 2 (~98 GB), `batch_size=1024` (VRAM probe kept 1024; ~85 GB)

## Progress

| Item | Status |
|------|--------|
| Jobs | **1 / 30** — `stage0-12-nlight-oof-f0-r0` (seed 42) |
| Epoch | **13 / 16** (early stopping patience 5; val_loss not beating epoch 1) |
| Remaining | f0 r1–r5, then folds 1–4 |

Runner: `scripts/run_12_nlight_oof.sh` (launched `15:11` after dense-S1 queue
stopped). Dense S1 **fold 0** checkpoint kept; fold 1 was incomplete.

Nine-pack 3-fold nested reference (not this 5-fold OOF): **0.368 / 9.88 / 0.803**.
