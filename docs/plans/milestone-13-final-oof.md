# Milestone 13: Final study-grouped OOF (5×6)

> **Alias:** historical **Milestone 7** in ADR 0007 / older docs.
> Index: [`MILESTONE_INDEX.md`](MILESTONE_INDEX.md).

**Status:** `blocked` until Milestones **9–11** (gene-only selection, scale
campaign decisions, fold-selected full model) are complete enough to justify
the expensive protocol.

**Protocol:** 5 outer folds × up to 6 restarts; orientation-aligned scores
([ADR 0008](../adr/0008-score-identifiability.md)); no TBS
([ADR 0009](../adr/0009-drop-tbs-scores.md)).

Prerequisites narrative (still valid; numbers updated):
[ADR 0007](../adr/0007-crossfit-prerequisites.md),
[`post-v0-scientific-programme.md`](post-v0-scientific-programme.md).

A 3-fold / 1-restart smoke of existing machinery is allowed for plumbing only
and must not overwrite v0.1 freezes.
