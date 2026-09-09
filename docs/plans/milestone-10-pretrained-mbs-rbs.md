# Milestone 10: Pretrained MBS/RBS scale campaign

> **Alias:** historical **7H**. Index: [`MILESTONE_INDEX.md`](MILESTONE_INDEX.md).

**Status:** `in_progress` (GPU 2 N-light OOF; GPU 0 on 12b probes)

Authoritative running log:
[`milestone-7h-pretrained-mbs-rbs-campaign.md`](milestone-7h-pretrained-mbs-rbs-campaign.md).
Latest architecture-decision detail:
[`reports/inspection/stage0_7h_nine_pack_smoke/vector_vs_scalar.md`](../../reports/inspection/stage0_7h_nine_pack_smoke/vector_vs_scalar.md).

## Sub-tracks

| ID | Status | Intent |
|----|--------|--------|
| **10a** | **done** | 5-combo pooling grid complete; **P2-G scalar max/max locked as cascade finalist** (no longer provisional) |
| **10a-cap** | **done (adopted)** | N-light `rho_hidden` 10→64 **adopted as default** (`mbs_e2e` 0.308/14.7/0.852); stale `external_test` 0.194/20.6 was a wrong readout |
| **10a-age-cov** | done (rejected) | Age head tissue+sex conditioning — full 3-fold ablation, all 3 metrics worse, not adopted |
| **10a-warmstart** | **done** | Vector (`region_hidden`) warm-started from P2-G/scalar-mean-max — both max/max (0.342/13.41/0.851) and mean/max (0.345/14.58/0.862) closed most/all of the from-scratch gap; scalar stays finalist |
| **10a-rbs-vs-mbs** | **resolved → 10e** | `rbs_linear_probe` beats `mbs_e2e`; OOF co-primary is frozen **enet** on RBS/MBS, not joint e2e alone |
| **10a-dense-stage1** | **partial** | Fold 0 done (e2e 0.348/14.52, linear age 10.93); queue stopped mid fold 1 to free GPU 2 for N-light OOF |
| **10a-onehop-smoke** | **done, both pass** | One-hop correctness smokes — first-ever full completion (3rd bug found+fixed). G0 beats G1 (0.335/17.3/0.750 vs 0.146/21.7/0.560), matching 9c. Multi-seed 42/43/44 genuinely distinct (F1 span 0.027, not collapsed) — restart ensembling is trustworthy for Milestone 12 |
| **10c** | **done** | Cancer pack **0.954**; Alzheimer’s **0.838**; subtype macro-F1 **0.202**; BMI not useful |
| **10d** | pending | Reference checkpoints after native P2-G + enet readout |
| **10e** | **done (negative)** | S1–S4 1-fold loses to native P2-G; staged recipe not promoted — see [`milestone-10e-staged-rbs-mbs-training.md`](milestone-10e-staged-rbs-mbs-training.md) |
| deferred | post-recipe | Sex-chromosome sex imputation, epigenetic-clock age imputation — after Milestone 12 / 12b recipe |
| **12b probe** | **queued GPU0** | CpGPT 65k N-light ablation (`stage0_12b_cpgpt_nlight_smoke`) — isolated probe after dense gene-expansion smoke; **not** “strictly after Milestone 12” and **not** the 12b product recipe |

## Priority rationale

- **P2-G scalar max/max is now the locked cascade finalist** — full 5-combo
  pooling grid at nine-pack scale settled it (see `vector_vs_scalar.md`).
- **Freeze-and-reuse is the standing architecture philosophy going forward**
  (DeepRVAT precedent: train a gene-invariant encoder once, freeze it, reuse
  frozen scores for new traits via lightweight probes rather than joint
  retraining). Applies to: BMI/ancestry and diagnoses **below** the n>200
  encoder-aux bar; n>200 disease/cancer classes now supervise the OOF
  encoder. Also applies to warm-starting the vector aggregation variant
  from the converged scalar encoder.
- **BMI/ancestry labels are joined** (2 070 / 1 380); heads stubbed for
  freeze-reuse. Ancestry has no class ≥1k. **Brain/blood** still have no
  honest trait head (catalogue / control-only). Disease classes and cancer
  types with **n>200** nine-pack disease-tissue samples are **encoder aux
  heads** on Milestone 12 (not freeze-reuse-only). Remaining diagnoses stay
  freeze-reuse / census.
- **Milestone 12 N-light** is running (30-ep + n>200 aux). **10e staged
  recipe rejected** — cascade 5×6, if launched, is **native P2-G** with
  frozen **enet** readouts. Extra Hub packs / 173k samples = freeze-reuse
  traits, not a bigger joint train.

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

Do not auto-launch Milestone **11** Stage B GPU or **cascade** 5×6 from this
campaign. N-light 5×6 is started from `scripts/run_12_nlight_oof_gpu2_after_s1.sh`.
This OOF cohort is HM450; the encoder must stay gene-invariant for later
EPIC / ONT (more probes per gene) — no mixed-platform claim yet.
