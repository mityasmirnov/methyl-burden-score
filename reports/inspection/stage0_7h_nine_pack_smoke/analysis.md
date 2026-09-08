# Milestone 10 — nine-pack campaign analysis

Updated: `2026-09-08T14:16+02:00`

- Matrix: `matrix-hub-nine-pack-virtual-v1` · Split: `hub-nine-pack-3fold-v1` (**34 234**)
- Platform: **HM450 only** · Primary: **`mbs_e2e`** outer test
- Detail dump: [`vector_vs_scalar.md`](vector_vs_scalar.md)

---

## Live status

| Job | Status | Detail |
|-----|--------|--------|
| Vector warm-starts (max + mean) | **done** | See leaderboard |
| **One-hop smokes** | **running** GPU 0 | G0 ≈0.335 / 17.3 / 0.75; G1 in flight |
| Dense stage-1 (optional) | other GPU | not a finalist blocker |
| Milestone 11 Stage B | deferred parallel | does **not** block OOF |

---

## Executive decisions

| Decision | Status |
|----------|--------|
| Cascade finalist | **P2-G scalar max/max** (0.355 / 13.431 / 0.853) |
| Light finalist | **N-light@64** (0.308 / 14.727 / 0.852) |
| Cold / warm vector | Warm closes gap; **neither replaces P2-G** |
| Age covariates | **Rejected** |
| Milestone **12** | OOF = those two finalists; needs one-hop smokes; **not gated on 11** |

---

## Leaderboard (`mbs_e2e`)

| Arm | Tissue F1 | Age MAE | Sex AUROC |
|-----|----------:|--------:|----------:|
| **P2-G scalar max/max** | **0.355** | **13.431** | 0.853 |
| vector mean→max **warm** | 0.345 | 14.578 | 0.862 |
| vector max→max **warm** | 0.342 | **13.406** | 0.851 |
| scalar mean/max | 0.330 | 13.496 | **0.875** |
| vector mean→max cold | 0.335 | 16.448 | 0.780 |
| vector max→max cold | 0.333 | 15.408 | 0.834 |
| scalar max/mean | 0.318 | 17.426 | 0.812 |
| **N-light@64** | 0.308 | 14.727 | 0.852 |
| N-light@10 | 0.273 | 15.057 | 0.742 |
| P2-G + age covariates | 0.340 | 14.135 | 0.834 |

**Interpretation.** P2-G wins the cold grid on the weighted tasks. Warm vector ≈ parity on age for max→max but still −0.013 tissue; mean→max warm has a bad fold-2 age. **Finalists unchanged.**

---

## Results digest (vector_vs_scalar)

1. **Cold 5-combo:** P2-G best tissue+age; cold vector loses.  
2. **rho 10→64:** adopted for N-light (+0.035 F1, −0.33 y, +0.11 AUROC).  
3. **Age covariates:** rejected.  
4. **Warm-starts:** both pooling recipes done; optimization helps; no finalist flip.  
5. **Side note:** `rbs_linear_probe` often beats `mbs_e2e` on age — discuss for OOF reporting, don’t change selection unilaterally.

---

## One-hop smokes (in flight)

G0 (dense, no mask): tissue **0.335**, age **17.3**, sex **0.750** (5 ep ATS fold 0, rho=64).  
G1 (seed-mask) + multi-seed 42/43/44 next.  
Report: `reports/inspection/stage0_7h_onehop_correctness/`.

**Interpretation.** Cheap correctness before Milestone 12 light-model assumptions.

---

## Outlook / next steps

1. Finish one-hop G1 + multi-seed → publish smoke report.  
2. Schedule **Milestone 12** 5×6 OOF on **P2-G + N-light@64** (gene-linked).  
3. Keep **Milestone 11** as parallel sparse-panel work (optional).  
4. Optional: dense stage-1 / RBS co-readout — science follow-ups only.  
5. Soft-stop: no full-arm OOF, no pack-mask trait GPU, no auto Stage B from keeper.

---

## Trait hygiene (10c)

age/sex/tissue wired · disease/cancer `label_status` ready · BMI/ancestry joined/stubbed · blood/brain defer.  
See [`trait_hygiene.md`](trait_hygiene.md).
