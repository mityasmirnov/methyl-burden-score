# Milestone 12: Final study-grouped OOF (5×6)

> **Alias:** historical **Milestone 7** in ADR 0007 / older docs.
> Index: [`MILESTONE_INDEX.md`](MILESTONE_INDEX.md).
>
> **Renumber note (2026-09-07):** previously stubbed as Milestone **13**;
> expression auxiliary is now **13** (after OOF).

**Status:** `blocked` until Milestone **10** is sufficient (**+ one-hop
correctness smokes**). **Does not wait on Milestone 11.**

**Protocol:** 5 outer folds × up to 6 restarts; orientation-aligned scores
([ADR 0008](../adr/0008-score-identifiability.md)); no TBS
([ADR 0009](../adr/0009-drop-tbs-scores.md)).

## Arms (finalists only)

After Milestone **10** scale screen:

1. **Cascade** — **P2-G scalar max/max** (locked). Vector warm-start improved
   cold vector (max→max warm 0.342 / 13.406 / 0.851) but did not beat P2-G
   tissue (0.355); not a finalist change unless mean→max warm clearly wins.
2. **Light** — **N-light@64** one-hop (`m_only` / gene-mean) — co-equal
   deployment model (`rho_hidden=64` default).

Do **not** spend 5×6 budget on all ~9 Stage A arms.

## Panel / cohort

- **This OOF:** gene-linked path used in Milestone **10** nine-pack / ATS refs
  (same setting that selected the finalists).
- **Not required:** Milestone **11** fold-selected sparse panel. A separate
  sparse-panel OOF may be scoped later after Stage B; it is not a prerequisite.

## Prerequisites

| Required | Optional / parallel |
|----------|---------------------|
| Milestone **9** gene-only selection | Milestone **11** Stage B (sparsity product) |
| Milestone **10** scale + finalists (P2-G, N-light@64) | Freeze-reuse trait probes (10c) |
| One-hop correctness smokes (seed-mask + multi-seed) | |

Narrative: [ADR 0007](../adr/0007-crossfit-prerequisites.md),
[`post-v0-scientific-programme.md`](post-v0-scientific-programme.md).
Campaign numbers:
[`../reports/inspection/stage0_7h_nine_pack_smoke/analysis.md`](../reports/inspection/stage0_7h_nine_pack_smoke/analysis.md).

A 3-fold / 1-restart smoke of existing machinery is allowed for plumbing only
and must not overwrite v0.1 freezes.

**Next:** Milestone **13** — expression continue-training / finetune (optional).
