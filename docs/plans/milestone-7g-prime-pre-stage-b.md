# Plan: 7G′ pre–Stage-B reporting + typed RBS + seed-gene Stage A gate

Status: **reporting + R0–R5 done** — **no P2-G architecture lock**. Current
GPU gate is the **matched 16-epoch promotion screen**
([`milestone-7g-prime-16ep-promotion.md`](milestone-7g-prime-16ep-promotion.md)).
Age-primary seed-mask scaffolding + fold-0 audit are **done**; CUDA waits on
16-ep unlock
([`milestone-7g-prime-age-seed-mask.md`](milestone-7g-prime-age-seed-mask.md)).
Fold-selected-panel Stage B is later and separate. R0–R5 does **not** decide
Stage B. Neural typed pool **not** promoted (R1–R3 age MAE ↑ vs R0, but
within-gene role shuffle did not collapse).

Parent: [`milestone-7g-prime-stage-a-deeprvat-screen.md`](milestone-7g-prime-stage-a-deeprvat-screen.md).
Report: [`reports/inspection/stage0_7g_gene_only_probe/analysis.md`](../../reports/inspection/stage0_7g_gene_only_probe/analysis.md).

## Scope and acceptance

1. Trustworthy Stage A reporting (no outer-test peeks; metric-specific fold
   counts; actual best epoch / samples seen / optimizer steps; graph content
   hash; honest N-light epoch labels).
2. Drop any declared/retained P2-G (or cascade) architecture lock from the ATS
   screen.
3. CPU typed-RBS ablation R0–R5 on saved scores (concurrent, cheap).
4. Retarget the next gate to a trait/seed-gene **Stage A repeat**.

**Done when:** report + `lock_recommendation.json` show `architecture_locked:
false`; analysis states the seed-gene gate; R0–R5 summary exists under
`reports/inspection/stage0_7g_gene_only_probe/typed_rbs_pooling/`; unit tests
cover fold counts and typed pooling.

## Locked decisions

| Choice | Decision | Why |
|--------|----------|-----|
| Architecture lock | **None** from ATS screen | Tissue-primary all-gene screen ≠ DeepRVAT seed protocol |
| Next gate | Trait/seed-gene Stage A repeat | Real DeepRVAT-analogue selection |
| Stage B fold-panel | Later, separate | Not unblocked by R0–R5 or free GPU |
| R0–R5 | Concurrent CPU diagnostic | Cheap; pooling evidence only |
| Neural typed pool | Follow-up only if R1–R4 beat R0 **and** shuffle collapses | Not a Stage B go/no-go; **blocked** after R0–R5 (shuffle held) |

## Seed-gene Stage A repeat (next gate)

Relabel the former “Stage B seed-gene transfer / design only” sketch as a
**Stage A repeat**:

1. Join EWAS Atlas CpG→trait tables to Hub sample/study IDs (start from
   `reports/inspection/deepmat_data_v1/trait_eligibility.md`; do not invent a trait).
2. Pick **one** trait that clears the bar (labels, Atlas CpG set, genes after
   `explicit_only`, ≥2 independent studies). If none, stop and document.
3. Persist seed gene IDs with study-overlap control.
4. Re-run the Stage A architecture screen with the shared encoder trained on
   **seed genes only**, applied to all eligible genes; evaluate non-seed transfer.

## Non-goals

- Declaring a P2-G lock
- Starting fold-panel Stage B because R0–R5 finished or GPU 0 is free
- 30-epoch N-light retrain / gated one-hop training in this change
- Neural typed aggregator in this change
- Raising LR
