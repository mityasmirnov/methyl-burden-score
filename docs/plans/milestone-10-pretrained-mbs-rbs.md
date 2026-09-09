# Milestone 10: Pretrained MBS/RBS scale campaign

> **Alias:** historical **7H**. Index: [`MILESTONE_INDEX.md`](MILESTONE_INDEX.md).

**Status:** `in_progress` (GPU 2 N-light OOF ← NOW; GPU 0 on 12b probes;
cascade blocked on GATE G1–G4)

Authoritative running log:
[`milestone-7h-pretrained-mbs-rbs-campaign.md`](milestone-7h-pretrained-mbs-rbs-campaign.md).
Latest architecture-decision detail:
[`reports/inspection/stage0_7h_nine_pack_smoke/vector_vs_scalar.md`](../../reports/inspection/stage0_7h_nine_pack_smoke/vector_vs_scalar.md).

## Sub-tracks

| ID | Status | Intent |
|----|--------|--------|
| **10a** | **done** | 5-combo pooling grid complete; **P2-G scalar max/max locked as cascade *topology* finalist** |
| **10a-cap** | **done (adopted)** | N-light `rho_hidden` 10→64 **adopted as default** (`mbs_e2e` 0.308/14.7/0.852); stale `external_test` 0.194/20.6 was a wrong readout |
| **10a-age-cov** | done (rejected) | Age head tissue+sex conditioning — full 3-fold ablation, all 3 metrics worse, not adopted |
| **10a-warmstart** | **done** | Vector (`region_hidden`) warm-started from P2-G/scalar-mean-max — both max/max (0.342/13.41/0.851) and mean/max (0.345/14.58/0.862) closed most/all of the from-scratch gap; scalar stays finalist |
| **10a-rbs-vs-mbs** | **resolved → 10e** | `rbs_linear_probe` beats `mbs_e2e`; OOF co-primary is frozen **enet** on RBS/MBS, not joint e2e alone |
| **10a-dense-stage1** | **partial** | Fold 0 done (e2e 0.348/14.52, linear age 10.93); queue stopped mid fold 1 to free GPU 2 for N-light OOF |
| **10a-onehop-smoke** | **done, both pass** | One-hop correctness smokes — G0 beats G1; multi-seed distinct — restart ensembling trustworthy for Milestone 12 |
| **10c** | **done** | Cancer pack **0.954**; Alzheimer’s **0.838**; subtype macro-F1 **0.202**; BMI not useful |
| **10d** | pending | Reference checkpoints — **GATE G4** after fair 10e + OOF finalist |
| **10e** | **in_progress / reopen** | Truncated S1→LP-FT ≠ full recipe; fair S1→S2→S3→S4 required — **GATE G4** |
| deferred | post-GATE | Sex-chromosome sex imputation, epigenetic-clock age imputation |
| **G2 probe** | **queued GPU0** | CpGPT 65k N-light ablation — GATE G2; not the 12b product recipe |

## Priority rationale

- **P2-G scalar max/max is the locked cascade *topology* finalist** — nine-pack
  5-combo grid (see `vector_vs_scalar.md`).
- **Freeze-and-reuse** remains the standing philosophy for traits below the
  n>200 encoder-aux bar; n>200 disease/cancer classes supervise Milestone 12
  encoder aux heads.
- **Milestone 12 N-light** is NOW. **Cascade 5×6 is blocked on GATE G1–G4**
  (12b gene util, CpGPT/positional, 12c platform, fair 10e + 10d). Default
  cascade recipe stays **native P2-G** until fair 10e flips it. Extra Hub
  packs / 173k samples = freeze-reuse traits, not a bigger joint train.

## Runners

- Warm-start + age-covariates queue: `scripts/run_7h_age_cov_and_warmstart_queue.sh` (done)
- Dense stage-1 pretrain + transplant queue (GPU 2): `scripts/run_7h_dense_stage1_queue.sh`
- One-hop correctness smokes (GPU 0): `scripts/run_7h_onehop_correctness_smokes.py`
- Vector configs: `configs/experiment/stage0_7h_nine_pack_vector_*.yaml`
  (`*_warmstart.yaml` = LP-FT from P2-G/scalar-mean-max; `*_from_dense_stage1.yaml`
  = LP-FT from the dense-gradient stage-1 checkpoint)
- Age-covariates config: `configs/experiment/stage0_7h_nine_pack_p2_g_age_covariates.yaml`
- Label status: `scripts/write_nine_pack_label_status.py`

## Hard stops

Do not auto-launch **cascade** 5×6 before GATE G1–G4 close (or are waived in
`TODO_PIPELINE.md`). Do not auto-launch full Milestone **11** Stage B GPU from
this campaign. N-light 5×6 owns GPU 2
(`scripts/run_12_nlight_oof_gpu2_after_s1.sh`). This OOF cohort is HM450; the
encoder must stay gene-invariant for later EPIC / ONT — no mixed-platform
claim yet (see **12c / G3**).
