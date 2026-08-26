# Stage 0 Milestone 7F — RBS→gene + direct leftover

Status: **done** for topology acceptance (assignment + trainer + saved-score
fusion + inspection). Full-budget methylation-only bake-off is **7G**.

## Topology (ADR 0009)

```text
CpG → typed region (gene roles | RBS) → RBS score
        ├─ allocated to gene (typed and/or nearest-gene) → MBS
        └─ no gene allocation → orphan RBS
CpG with no typed region (incl. former TBS-only) → direct
late fusion: [orphan RBS | MBS | direct] → linear heads
```

- **No TBS arm** in the model matrix or score export.
- Nearest-gene applies only to typed **RBS** allocation, never leftover CpGs.

## Artifacts

| Item | Path |
|------|------|
| Plan | `docs/plans/milestone-7f-rbs-gene-direct.md` |
| ADR | `docs/adr/0009-drop-tbs-scores.md` |
| Config | `configs/experiment/stage0_7f_rbs_gene_direct.yaml` |
| CLI | `mbs train cascade` |
| Fixture scores | `artifacts/runs/stage0-7f-fixture-v1/scores/` |
| Hub smoke (fold 0) | `artifacts/runs/stage0-7f-hub-smoke-v1/fold_0/scores/` |
| Frozen folds | `hub-ats-7e-3fold-v1` |

Score dirs contain `mbs.zarr`, `rbs.zarr` (orphan), `direct_contrib.zarr`,
`score_manifest.json` — **no** `tbs.zarr`.

## Evidence

1. **Unit fixtures** (`tests/unit/test_stage0_7f.py`): leftover→direct;
   RBS→gene nearest-gene; TBS edges ignored; fusion rejects TBS keys.
2. **Fixture end-to-end**: `mbs train cascade --overfit-fixture` wrote this
   report and fused saved score matrices.
3. **Hub smoke on frozen folds**: one outer fold of `hub-ats-7e-3fold-v1`
   with documented smoke budget (`max_loci=512`, `max_epochs=1`,
   `max_train_samples=48`). See `summary.json` (`tbs_arm: false`).

Default Hub config retains the 7E ceiling (8192 loci / 2 epochs); raise budget
in **7G**.
