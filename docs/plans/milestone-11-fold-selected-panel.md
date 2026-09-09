# Milestone 11: Fold-selected panel + full model

> **Alias:** historical **7G′ Stage B**. Index: [`MILESTONE_INDEX.md`](MILESTONE_INDEX.md).

**Status:** `deferred` (parallel / optional — **does not hard-block Milestone 12**)

Parent brief:
[`milestone-7g-prime-matched-probe-lightweight.md`](milestone-7g-prime-matched-probe-lightweight.md)
(Stage B sections). Pre-gate notes:
[`milestone-7g-prime-pre-stage-b.md`](milestone-7g-prime-pre-stage-b.md).
Related: [`milestone-7h-fold-safe-probe-panel-benchmark.md`](milestone-7h-fold-safe-probe-panel-benchmark.md).

**Runner:** `scripts/run_7g_prime_stage_b.py`

## Policy (2026-09-08)

Milestone **10** locked the architecture finalists (**P2-G** cascade +
**N-light@64**). Milestone **11** is no longer the architecture gate and **no
longer blocks** Milestone **12** finalist OOF on the gene-linked path.

**11 remains valuable** as a **sparsity / fold-safe panel product** track:
classical `C-mvalue-enetS` on outer-train-selected loci, plus finalists retrained
on that same panel (`N-cascade-S`, light-on-S), optional fusion ablations, and
`direct_cpg.zarr` when direct loci exist.

- CPU `--panels-only` / `--classical-only` may run anytime.
- Full Stage B GPU: schedule explicitly after review; prefer a **thinner** arm
  matrix (panels + classical + finalists-on-S), not the full Stage A soup.
- Do **not** auto-launch Stage B GPU from the Milestone 10 keeper.

## Done when

- Fold-safe panels under
  `reports/inspection/stage0_7g_prime_matched_probe/fold_panels/fold_*_panel.json`
- Matched `C-mvalue-enetS` / finalist-on-S (`N-cascade-S`, light) / optional fusion
- `direct_cpg.zarr` when direct loci exist (`n_direct > 0`)
- Report under `reports/inspection/stage0_7g_prime_matched_probe/`

**Does not block:** Milestone **12** (gene-linked finalist OOF).

**12b gene-set option:** a fold-selected gene list (this milestone / ADR 0012:
discovery CpGs → seed genes → all linked CpGs of those genes) is a valid
DeepRVAT-style training-set alternative to “all ~19.6k represented genes” in
[`milestone-12b-full-gene-panel.md`](milestone-12b-full-gene-panel.md) §2c.
Still deferred; still not a gate on the 65k Milestone **12** 5×6.

## Arm matrix (identical panel per fold)

From [`configs/experiment/stage0_7g_prime_stage_b.yaml`](../../configs/experiment/stage0_7g_prime_stage_b.yaml):

| Arm | Role | Priority |
|-----|------|----------|
| `C-mvalue-enetS` | Sparse linear on fold-selected panel | **core** |
| `N-cascade-S` | P2-G cascade params on same panel | **core** (finalist-on-S) |
| `N-light-type` / gene-mean | One-hop on same panel | **core** (finalist-on-S) |
| `N-mbs-posthoc-full-fusion` | Orphan RBS + MBS + direct late fusion | optional |
| `N-mbs-posthoc-mbs-direct` | MBS + direct late fusion (no orphan) | optional |

Panel: `max_seeds: 10000`, split `hub-ats-7e-3fold-v1`, matrix
`matrix-hub-age-tissue-sex-full-v1`, graph-v2.

## CPU prep flags (no GPU)

```bash
# Write fold_panels only (CPU dry-run may use --panel-repeats 1; production: 5)
uv run python -u scripts/run_7g_prime_stage_b.py --device cpu --panels-only --folds 0 --panel-repeats 1

# Classical enetS after panels exist
uv run python -u scripts/run_7g_prime_stage_b.py --device cpu --classical-only --folds 0

# Or: wait for fold-0 panel, then classical + folds 1–2 (BLAS threads capped)
PANEL_REPEATS=1 bash scripts/run_stage_b_cpu_chain.sh
```

Cap OpenBLAS when selecting on wide matrices:
`OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8` (set in the chain script).
Stability selection univariate-prefilters to **4096** columns before the enet
grid (same screen as gene-seed panels); meta records `n_cols_after_prefilter`.

`--panels-only` and `--classical-only` are mutually exclusive. Panels-only
writes `fold_panels/manifest.json` and does **not** regenerate the full Stage B
report. Classical-only loads existing `fold_*_panel.json` and writes
`per_arm/C-mvalue-enetS.json`. `--panel-repeats` sets stability-selection
repeats (default 5; dry-run `1` is recorded in the panel JSON/manifest).
Lock input:
`reports/inspection/stage0_7g_gene_only_probe/lock_recommendation.json`
(update when promoting: cascade finalist **P2-G**, light **N-light@64**,
`next_gate` no longer “must finish before OOF”).

## Product export

When `assignment.n_direct > 0` and `locus_ids` are passed, cascade score writes
include `direct_cpg.zarr` + `direct_locus_index.parquet` (association product).
Gene-only (`assignment_gene_linked_only`) zeros direct columns — expected omit.
Orphan fusion needs **qualified multi-CpG** orphans; census:
`reports/inspection/stage0_7g_prime_matched_probe/orphan_rbs_census.md`.

## Non-goals

- Re-opening Milestone **10** architecture selection.
- Blocking Milestone **12** gene-linked finalist OOF.
- Auto-launch from the nine-pack GPU keeper.
