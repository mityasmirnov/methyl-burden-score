# Milestone 12: Final study-grouped OOF (5×6)

> **Alias:** historical **Milestone 7** in ADR 0007 / older docs.
> Index: [`MILESTONE_INDEX.md`](MILESTONE_INDEX.md).
>
> **Renumber note (2026-09-07):** previously stubbed as Milestone **13**;
> expression auxiliary is now **13** (after OOF).

**Status:** `blocked` until Milestones **9–11** (gene-only selection, scale
campaign decisions, fold-selected full model) are complete enough to justify
the expensive protocol.

**Protocol:** 5 outer folds × up to 6 restarts; orientation-aligned scores
([ADR 0008](../adr/0008-score-identifiability.md)); no TBS
([ADR 0009](../adr/0009-drop-tbs-scores.md)).

**Arms:** **finalists only** after Milestone **10** scale screen:

1. Winning **cascade** (scalar P2-G *or* nine-pack vector RBS — undecided until
   `reports/inspection/stage0_7h_nine_pack_smoke/vector_vs_scalar.md`)
2. **N-light** one-hop light MBS (`m_only` / gene-mean) — co-equal deployment
   model, not a side ablation

Do **not** spend 5×6 budget on all ~9 Stage A arms.

Prerequisites narrative (still valid; numbers updated):
[ADR 0007](../adr/0007-crossfit-prerequisites.md),
[`post-v0-scientific-programme.md`](post-v0-scientific-programme.md).

A 3-fold / 1-restart smoke of existing machinery is allowed for plumbing only
and must not overwrite v0.1 freezes.

**Next:** Milestone **13** — expression continue-training / finetune (optional).
