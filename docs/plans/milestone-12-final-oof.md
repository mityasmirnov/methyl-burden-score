# Milestone 12: Final study-grouped OOF (5×6)

> **Alias:** historical **Milestone 7** in ADR 0007 / older docs.
> Index: [`MILESTONE_INDEX.md`](MILESTONE_INDEX.md).
>
> **Renumber note (2026-09-07):** previously stubbed as Milestone **13**;
> expression auxiliary is now **13** (after OOF).

**Status:** `blocked` until Milestone **10e** staged RBS→MBS recipe is
smoked (**+ one-hop smokes + enet co-primary**). **Does not wait on
Milestone 11.** **Do not** launch 5×6 on today’s joint `mbs_e2e` P2-G /
N-light pair.

Recipe: [`milestone-10e-staged-rbs-mbs-training.md`](milestone-10e-staged-rbs-mbs-training.md).

**Protocol:** 5 outer folds × up to 6 restarts; orientation-aligned scores
([ADR 0008](../adr/0008-score-identifiability.md)); no TBS
([ADR 0009](../adr/0009-drop-tbs-scores.md)).

## Arms (finalists only) — recipe, not just architecture names

Architecture *candidates* after Milestone **10** still P2-G cascade topology
+ N-light@64, but **training and readout change before OOF**:

1. **Cascade** — staged **S1–S4** (dense vector RBS → freeze → gene MBS hop →
   unfreeze FT), not 15-ep joint e2e max/max. Report **`rbs_enet` +
   `mbs_enet` co-primary**; `mbs_e2e` secondary.
2. **Light** — N-light@64 encoder + **`mbs_enet_nested` as product readout**
   (ATS nested enet age 10 vs e2e 17). Do not 5×6 e2e-only.

Do **not** spend 5×6 budget on all ~9 Stage A arms or on the old joint-e2e
recipe. Extra traits (BMI, ancestry, cancer types, AD/PD/stroke, …) are
**frozen-score probes**, not extra OOF arms.

## Panel / cohort

- **This OOF:** gene-linked path used in Milestone **10** nine-pack / ATS refs
  (same setting that selected the finalists).
- **Not required:** Milestone **11** fold-selected sparse panel. A separate
  sparse-panel OOF may be scoped later after Stage B; it is not a prerequisite.

## Prerequisites

| Required | Optional / parallel |
|----------|---------------------|
| Milestone **9** gene-only selection | Milestone **11** Stage B (sparsity product) |
| Milestone **10** topology (P2-G, N-light@64) | Freeze-reuse trait probes (10c) — not extra OOF arms |
| **10e** 1-fold staged-recipe smoke + enet co-primary | |
| One-hop correctness smokes (seed-mask + multi-seed) | |

Narrative: [ADR 0007](../adr/0007-crossfit-prerequisites.md),
[`post-v0-scientific-programme.md`](post-v0-scientific-programme.md).
Campaign numbers:
[`../reports/inspection/stage0_7h_nine_pack_smoke/analysis.md`](../reports/inspection/stage0_7h_nine_pack_smoke/analysis.md).

A 3-fold / 1-restart smoke of existing machinery is allowed for plumbing only
and must not overwrite v0.1 freezes.

**Next:** Milestone **13** — expression continue-training / finetune (optional).
