# 12c. Platform robustness — GATE G3

> **Status:** `pending` (post–N-light **GATE G3**; blocks product cascade).
> Parent checklist: [`../TODO_PIPELINE.md`](../TODO_PIPELINE.md).
> Index: [`MILESTONE_INDEX.md`](MILESTONE_INDEX.md).
> Related: [`milestone-12b-full-gene-panel.md`](milestone-12b-full-gene-panel.md) §2b
> (within-gene invariance; no HM450 identity crosswalk).

## Question

Can training stay **gene-invariant** under platform-like missingness (random
CpG drop, 450K-like vs EPIC-like masks, varying counts per gene), and what is
the minimal path to consume EPIC (and later ONT) without treating absent
CpGs as low burden?

## Approaches tested

1. **Prefer first (no EPIC matrix required):** HM450 training-time CpG /
   platform-mask dropout — simulate EPIC↔450K-style coverage on the nine-pack
   matrix; evaluate MBS stability after downsampling.
2. **Follow-on when data ready:** build per-probe gene/region **membership**
   for EPIC’s own coordinates/manifest (graph join *or* Illumina gene columns);
   optional train mix. **Not** an HM450↔EPIC probe-ID crosswalk.
3. **Later:** ONT measurement adapter (frequency/depth/uncertainty) — Milestone
   14f / 12b follow-on D; not G3 acceptance.

## Results

`pending` — not started.

## Verdict

Open. Waive only with an explicit written verdict in `TODO_PIPELINE.md`.

## Done when

- Dropout / mask smoke report under `reports/inspection/` with nested or e2e
  comparison vs unmatched baseline.
- Short invariance note: membership label + observed-only pooling; unobserved
  genes stay `mbs_present=0`.
- EPIC membership table (if attempted) documented with source columns / graph id.

## Non-goals

- Mixing EPIC/ONT into the in-flight 65k N-light 5×6.
- Implementing the DeepRVAT sampler (that is **12b / G1**).
- Claiming platform-agnostic product without G3 evidence.
