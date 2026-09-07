# Plan: 7G′ pre–Stage-B reporting + typed RBS + seed-gene Stage A gate

> **Canonical milestones: 9d** (typed-RBS) / gate notes for **11**. Historical 7G′ pre–Stage-B. See [`MILESTONE_INDEX.md`](MILESTONE_INDEX.md).


Status: **reporting + R0–R5 done** — **no P2-G architecture lock**. Seed-mask
Stage A (9c) is **done** and **not adopted**. Matched 16-ep promotion (9b) is
**done**. Current scientific gate is Milestone **10** scale review
(`next_gate: milestone_10_scale_review` in
`lock_recommendation.json`). Fold-selected-panel Stage B is **Milestone 11**
and stays blocked until 10a+10b are reviewed. R0–R5 does **not** decide
Stage B. Neural typed pool **not** promoted (R1–R3 age MAE ↑ vs R0, but
within-gene role shuffle did not collapse).

Parent: [`milestone-7g-prime-stage-a-deeprvat-screen.md`](milestone-7g-prime-stage-a-deeprvat-screen.md).
Report: [`reports/inspection/stage0_7g_gene_only_probe/analysis.md`](../../reports/inspection/stage0_7g_gene_only_probe/analysis.md).
Stage B plan: [`milestone-11-fold-selected-panel.md`](milestone-11-fold-selected-panel.md).

## Scope and acceptance

1. Trustworthy Stage A reporting (no outer-test peeks; metric-specific fold
   counts; actual best epoch / samples seen / optimizer steps; graph content
   hash; honest N-light epoch labels).
2. Drop any declared/retained P2-G (or cascade) architecture lock from the ATS
   screen.
3. CPU typed-RBS ablation R0–R5 on saved scores (concurrent, cheap).
4. ~~Retarget the next gate to a trait/seed-gene Stage A repeat.~~ **Superseded:**
   9c cleared the seed-mask path (not adopted); next gate is Milestone **10**.

**Done when:** report + `lock_recommendation.json` show `architecture_locked:
false`; analysis states the current next gate; R0–R5 summary exists under
`reports/inspection/stage0_7g_gene_only_probe/typed_rbs_pooling/`; unit tests
cover fold counts and typed pooling.

## Locked decisions

| Choice | Decision | Why |
|--------|----------|-----|
| Architecture lock | **None** from ATS screen | Tissue-primary all-gene screen ≠ DeepRVAT seed protocol |
| Provisional reference | **P2-G** (not locked) | Best landed cascade; YAML Stage B defaults while `locked_cascade_arm` is null |
| Next gate | Milestone **10** scale review | Seed-mask (9c) cleared; need nine-pack + 10b pooling honesty |
| Stage B fold-panel | Milestone **11**, later | Not unblocked by R0–R5 or free GPU |
| R0–R5 | Concurrent CPU diagnostic | Cheap; pooling evidence only |
| Neural typed pool | Follow-up only if R1–R4 beat R0 **and** shuffle collapses | Not a Stage B go/no-go; **blocked** after R0–R5 (shuffle held) |

## Seed-gene Stage A repeat (historical; gate cleared)

The former “Stage A seed-gene repeat” / age-primary seed-mask path ran as
Milestone **9c**. Result: **G0 beats G1–G3** — seed-masking **not adopted**.
Do not reopen unless a new hypothesis appears. See
[`milestone-7g-prime-age-seed-mask.md`](milestone-7g-prime-age-seed-mask.md).

## Non-goals

- Declaring a P2-G lock
- Starting fold-panel Stage B GPU because R0–R5 finished or GPU 0 is free
- 30-epoch N-light retrain / gated one-hop training in this change
- Neural typed aggregator in this change
- Raising LR
