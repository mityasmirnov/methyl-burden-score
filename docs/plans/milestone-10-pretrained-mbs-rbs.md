# Milestone 10: Pretrained MBS/RBS scale campaign

> **Alias:** historical **7H**. Index: [`MILESTONE_INDEX.md`](MILESTONE_INDEX.md).

**Status:** `in_progress` (GPU 0)

Authoritative running log:
[`milestone-7h-pretrained-mbs-rbs-campaign.md`](milestone-7h-pretrained-mbs-rbs-campaign.md).

## Sub-tracks

| ID | Status | Intent |
|----|--------|--------|
| **10a** | wrapping | Nine-pack P2-G done; m-only f0 full-budget repair |
| **10a+** | queued | Nine-pack **vector RBS** @15 ep (mean→max, max→max) vs provisional P2-G |
| **10a++** | queued | One-hop correctness: seed-mask G0/G1 + multi-seed (fold0, ~5 ep) |
| **10b** | queued | ATS seed-43 2×2 pooling (after 10a+/++); B.5 census done |
| **10c** | partial | `label_status` sidecar written; wire disease/cancer heads later |
| **10d** | pending | Reference checkpoints **after** finalists chosen |

## Priority rationale

- **P2-G is provisional** (`architecture_locked: false`). ATS scalar vs vector was
  within noise; nine-pack vector is the honest scale test.
- **Validate one-hop tricks cheaply** before any multi-trait / multi-restart
  light-model campaign: seed-mask on one-hop was never run (9c = cascade only);
  multi-seed on one-hop has zero prior successful runs after the seed-offset fix.
- **Do not OOF all 9 packs.** Usable now: age/tissue/sex. Disease/cancer need
  `label_status` (sidecar ready). Blood/brain deferred. BMI/ancestry not joined.
- **Milestone 12** = finalists only (winning cascade + N-light).
- **Demoted:** ATS light-mean s2 keep-busy (lower value than 10a+/++).

## Runners

- Keeper: `scripts/run_7h_next_queue.sh` (+ `run_7h_queue_handoff.sh`)
- Vector configs: `configs/experiment/stage0_7h_nine_pack_vector_*.yaml`
- One-hop smokes: `scripts/run_7h_onehop_correctness_smokes.py`
- Label status: `scripts/write_nine_pack_label_status.py`
- B.4 pooling: `scripts/run_7h_ats_pooling_s2.sh`

## Hard stops

Do not auto-launch Milestone **11** Stage B GPU or Milestone **12** OOF from
this campaign. HM450 only — no cross-platform claim yet.
