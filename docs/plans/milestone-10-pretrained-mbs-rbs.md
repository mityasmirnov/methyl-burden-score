# Milestone 10: Pretrained MBS/RBS scale campaign

> **Alias:** historical **7H**. Index: [`MILESTONE_INDEX.md`](MILESTONE_INDEX.md).

**Status:** `in_progress` (GPU 0)

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
| **10a-warmstart** | running | Vector (`region_hidden`) warm-started from P2-G checkpoint (LP-FT: freeze 4ep, fine-tune 11ep) — fold 0 promising (0.348/12.88/0.949), folds 1-2 in progress |
| **10c** | next | Freeze-and-reuse: disease/cancer frozen-feature logistic probes on the finalized encoder (labels already built; no GPU retrain needed) |
| **10d** | pending | Reference checkpoints **after** finalists chosen |
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
- **`brain`/`bmi`/`ancestry` packs have no real trait labels yet** (matrix
  routing only, no head exists) — freeze-and-reuse can't extend there
  without label-prep work first; out of scope for this campaign.
- **Milestone 12** = finalists only (winning cascade + N-light).

## Runners

- Warm-start + age-covariates queue: `scripts/run_7h_age_cov_and_warmstart_queue.sh`
- Vector configs: `configs/experiment/stage0_7h_nine_pack_vector_*.yaml`
  (`*_warmstart.yaml` = LP-FT from the P2-G/scalar-mean-max checkpoint)
- Age-covariates config: `configs/experiment/stage0_7h_nine_pack_p2_g_age_covariates.yaml`
- Label status: `scripts/write_nine_pack_label_status.py`

## Hard stops

Do not auto-launch Milestone **11** Stage B GPU or Milestone **12** OOF from
this campaign. HM450 only — no cross-platform claim yet.
