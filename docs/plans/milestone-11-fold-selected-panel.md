# Milestone 11: Fold-selected panel + full model

> **Alias:** historical **7G′ Stage B**. Index: [`MILESTONE_INDEX.md`](MILESTONE_INDEX.md).

**Status:** `blocked` (CPU prep in progress; GPU not launched)

Parent brief:
[`milestone-7g-prime-matched-probe-lightweight.md`](milestone-7g-prime-matched-probe-lightweight.md)
(Stage B sections). Pre-gate notes:
[`milestone-7g-prime-pre-stage-b.md`](milestone-7g-prime-pre-stage-b.md).
Related: [`milestone-7h-fold-safe-probe-panel-benchmark.md`](milestone-7h-fold-safe-probe-panel-benchmark.md).

**Runner:** `scripts/run_7g_prime_stage_b.py`

**Hard stop:** do **not** launch full Stage B GPU until Milestone **10a+10b**
results are reviewed. Seed-gene / seed-mask gate is **cleared** (9c done;
seed-mask **not adopted**). Remaining gate is honest **10** scale/architecture
refs.

## Done when

- Fold-safe panels under
  `reports/inspection/stage0_7g_prime_matched_probe/fold_panels/fold_*_panel.json`
- Matched `C-mvalue-enetS` / `N-cascade-S` / `N-light-type` / fusion ablations
- `direct_cpg.zarr` when direct loci exist (`n_direct > 0`)
- Report under `reports/inspection/stage0_7g_prime_matched_probe/`

**Blocks:** Milestone **12** final OOF.

## Arm matrix (identical panel per fold)

From [`configs/experiment/stage0_7g_prime_stage_b.yaml`](../../configs/experiment/stage0_7g_prime_stage_b.yaml):

| Arm | Role |
|-----|------|
| `C-mvalue-enetS` | Sparse linear on fold-selected panel |
| `N-cascade-S` | Cascade on same panel (provisional P2-G train params from YAML while lock arm is null) |
| `N-light-type` | Flat region / typed one-hop on same panel |
| `N-mbs-posthoc-full-fusion` | Orphan RBS + MBS + direct late fusion |
| `N-mbs-posthoc-mbs-direct` | MBS + direct late fusion (no orphan block) |

Panel: `max_seeds: 10000`, split `hub-ats-7e-3fold-v1`, matrix
`matrix-hub-age-tissue-sex-full-v1`, graph-v2.

## CPU prep flags (no GPU)

```bash
# Write fold_panels only
uv run python -u scripts/run_7g_prime_stage_b.py --device cpu --panels-only --folds 0

# Classical enetS after panels exist
uv run python -u scripts/run_7g_prime_stage_b.py --device cpu --classical-only --folds 0
```

`--panels-only` and `--classical-only` are mutually exclusive. Panels-only
writes `fold_panels/manifest.json` and does **not** regenerate the full Stage B
report. Classical-only loads existing `fold_*_panel.json` and writes
`per_arm/C-mvalue-enetS.json`.

Lock input:
`reports/inspection/stage0_7g_gene_only_probe/lock_recommendation.json`
(`architecture_locked: false`, `locked_cascade_arm: null`,
`provisional_reference_arm: P2-G`, `next_gate: milestone_10_scale_review`).

## Product export

When `assignment.n_direct > 0` and `locus_ids` are passed, cascade score writes
include `direct_cpg.zarr` + `direct_locus_index.parquet` (association product).
Gene-only (`assignment_gene_linked_only`) zeros direct columns — expected omit.
Orphan fusion needs **qualified multi-CpG** orphans; census:
`reports/inspection/stage0_7g_prime_matched_probe/orphan_rbs_census.md`.

## Non-goals

- Do not treat `_staging_N_cascade_S_fold_*` / `_staging_N_light_type_*` stubs as
  Stage B progress.
- Do not auto-launch from the Milestone 10 campaign queue.
- Do not rewrite phenotype tables here (10c policy is separate).
