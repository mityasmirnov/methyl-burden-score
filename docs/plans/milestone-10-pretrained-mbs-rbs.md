# Milestone 10: Pretrained MBS/RBS scale campaign

> **Alias:** historical **7H**. Index: [`MILESTONE_INDEX.md`](MILESTONE_INDEX.md).

**Status:** `in_progress` (GPU 0 + GPU 2 both loaded)

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
| **10a-dense-stage1** | running (GPU 2) | Dense-gradient (`cpg_pool`/`region_pool`: mean) stage-1 RBS pretrain → transplant into scalar + vector max/max, testing whether a more thoroughly-trained encoder beats P2-G's own result |
| **10a-onehop-smoke** | running (GPU 0) | One-hop correctness smokes (seed-mask + multi-seed) — 3rd bug found+fixed this relaunch, never completed before |
| **10c** | next | Freeze-and-reuse probes on frozen RBS/MBS (disease/cancer/BMI/ancestry/individual diagnoses) |
| **10d** | pending | Reference checkpoints **after** staged recipe + enet readout |
| **10e** | **pending (blocks 12)** | Staged RBS→MBS training + frozen enet co-primary — see [`milestone-10e-staged-rbs-mbs-training.md`](milestone-10e-staged-rbs-mbs-training.md) |
| deferred | post-OOF | CpGPT/positional embeddings, sex-chromosome sex imputation, epigenetic-clock age imputation — explicitly after Milestone 12 |

## Priority rationale

- **P2-G scalar max/max is now the locked cascade finalist** — full 5-combo
  pooling grid at nine-pack scale settled it (see `vector_vs_scalar.md`).
- **Freeze-and-reuse is the standing architecture philosophy going forward**
  (DeepRVAT precedent: train a gene-invariant encoder once, freeze it, reuse
  frozen scores for new traits via lightweight probes rather than joint
  retraining). Applies to: extending to disease/cancer (real labels exist,
  no head trained yet) and to warm-starting the vector aggregation variant
  from the converged scalar encoder.
- **BMI/ancestry labels are joined** (2 070 / 1 380); heads stubbed for
  freeze-reuse. Ancestry has no class ≥1k. **Brain/blood** still have no
  honest trait head (catalogue / control-only). Individual diseases
  (AD 945, PD 333, stroke 204, …) are freeze-reuse targets after a
  per-diagnosis n census — not joint encoder retrains.
- **Milestone 12** waits on **10e**: do not 5×6 joint `mbs_e2e` as-is.
  Staged RBS→MBS (dense vector RBS → freeze → MBS hop → unfreeze) plus
  frozen **enet** readouts (critical for N-light). Extra Hub packs / 173k
  samples = freeze-reuse traits, not a bigger joint train.

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

Do not auto-launch Milestone **11** Stage B GPU or Milestone **12** OOF from
this campaign. HM450 only — no cross-platform claim yet.
