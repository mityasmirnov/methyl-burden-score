# Milestone 12: Final study-grouped OOF (5×6)

> **Alias:** historical **Milestone 7** in ADR 0007 / older docs.
> Index: [`MILESTONE_INDEX.md`](MILESTONE_INDEX.md).
>
> **Renumber note (2026-09-07):** previously stubbed as Milestone **13**;
> expression auxiliary is now **13** (after OOF).

**Status:** **N-light arm `in_progress`** (GPU 2 after dense S1). **Cascade
arm `blocked`** until 10e S1–S4 smoke + P2-G nested-enet ranking.
**Does not wait on Milestone 11.**

Recipe / readout: [`milestone-10e-staged-rbs-mbs-training.md`](milestone-10e-staged-rbs-mbs-training.md).

**Protocol:** 5 outer folds × up to 6 restarts; orientation-aligned scores
([ADR 0008](../adr/0008-score-identifiability.md)); no TBS
([ADR 0009](../adr/0009-drop-tbs-scores.md)).

**Split (frozen):** `hub-nine-pack-5fold-v1` (34 234 samples, 470 studies).
Do **not** reuse `hub-nine-pack-3fold-v1` for this OOF.

**Product readout for every trait:** frozen nested elastic-net on encoder
scores (`mbs_enet_nested`; cascade also `rbs_enet`). Joint `mbs_e2e` is
encoder supervision / diagnostic only.

## Arms — N-light first

1. **Light (now)** — N-light@64 one-hop, mean pool, gene-linked CpGs.
   Runner: `scripts/run_12_nlight_oof.sh` (GPU 2, `batch_size` 1024 with
   VRAM probe). Run ids `stage0-12-nlight-oof-f{fold}-r{restart}`.
   Nine-pack 3-fold nested already **0.368 / 9.88 / 0.803** vs e2e
   **0.308 / 14.73 / 0.852** — do not 5×6 e2e-only.
2. **Cascade (later)** — P2-G topology trained **S1–S4**, not 15-ep joint
   max/max. Finalist vs N-light decided **under nested enet**, including
   whether vector `gene_rho` (`region_hidden`) beats scalar.

Extra traits with **n≥600** (BMI, ancestry, cancer pack, Alzheimer’s 945)
are **frozen-score probes**, not extra OOF arms. PD/stroke/etc. stay below
the case bar.

## Panel / cohort / platforms

- **This OOF:** gene-linked nine-pack HM450 (same encoder family as M10).
- **Architecture invariant:** ragged, gene-invariant Deep Set — must accept
  HM450, EPIC (more probes/gene), and ONT (even more) without a fixed
  manifest. This 5×6 does **not** mix platforms; transfer eval is separate.
- **Not required:** Milestone **11** fold-selected sparse panel.

## Prerequisites

| Required for N-light 5×6 | Required for cascade 5×6 | Optional / parallel |
|--------------------------|--------------------------|---------------------|
| Milestone **9** gene-only selection | 10e 1-fold S1–S4 smoke | Milestone **11** Stage B |
| N-light@64 + 3-fold nested enet | P2-G nested enet ranking | Freeze-reuse n≥600 traits |
| One-hop smokes (done) | | |
| `hub-nine-pack-5fold-v1` | | |

Narrative: [ADR 0007](../adr/0007-crossfit-prerequisites.md),
[`post-v0-scientific-programme.md`](post-v0-scientific-programme.md).

**Next:** Milestone **13** — expression continue-training / finetune (optional).
