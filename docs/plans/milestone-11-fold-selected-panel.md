# Milestone 11: Fold-selected panel + full model

> **Alias:** historical **7G′ Stage B**. Index: [`MILESTONE_INDEX.md`](MILESTONE_INDEX.md).

**Status:** `pending` for **GATE G1** gene-set decision (CPU prep anytime;
full Stage B GPU still scheduled explicitly). Does **not** hard-block
Milestone **12 N-light 65k**; **is in scope before product cascade**.

Parent brief:
[`milestone-7g-prime-matched-probe-lightweight.md`](milestone-7g-prime-matched-probe-lightweight.md)
(Stage B sections). Pre-gate notes:
[`milestone-7g-prime-pre-stage-b.md`](milestone-7g-prime-pre-stage-b.md).
Related: [`milestone-7h-fold-safe-probe-panel-benchmark.md`](milestone-7h-fold-safe-probe-panel-benchmark.md).

**Runner (legacy Stage B):** `scripts/run_7g_prime_stage_b.py`  
**Runner (scalable universe):** `scripts/run_trait_universe_catalog.py` +
`scripts/run_trait_universe_fold_reselect.py`  
**Config:** [`configs/experiment/stage0_11_trait_universe.yaml`](../../configs/experiment/stage0_11_trait_universe.yaml)

## Three-stage gene universe (2026-09-10)

The 65k elastic-net stability path does not scale to ~374k gene-linked CpGs ×
many traits. Replace discovery with **loose univariate EWAS** in three stages:

| Stage | Leakage (ADR 0011) | What |
|-------|--------------------|------|
| **A** catalog | `catalog_universe` | Gene-linked CpGs only; all-data univariate EWAS per trait; **union** EWAS Atlas genes |
| **B** restrict | `catalog_universe` | Keep CpGs of genes associated with **any** trait → working universe ≪ 482k |
| **C** fold-safe | `internal_fold` | Outer-train only: reselect + expand seeds → sibling gene CpGs |

**Invariant:** CV / models train on **C**, never on A. Folds reselect independently
(gene sets differ but may overlap). Encoder stays gene-permutation invariant.
Selection stays **loose** (nominal p &lt; 0.05 + gene floors). Disease/cancer:
case≠control via `label_status`; unknowns stay unknown (10c).

```bash
uv run python -u scripts/run_trait_universe_catalog.py \
  --config configs/experiment/stage0_11_trait_universe.yaml
uv run python -u scripts/run_trait_universe_fold_reselect.py \
  --config configs/experiment/stage0_11_trait_universe.yaml
```

Report: `reports/inspection/stage0_11_trait_universe/`.

## Question / Approaches / Results / Verdict

- **Question:** As a DeepRVAT-style gene panel, does fold-selected genes
  (ADR 0012) match or beat “all ~19.6k represented genes” under nested enet
  for product cascade?
- **Approaches tested:**
  - CPU Stage B: 3-fold ATS panels (`n_repeats=1` dry-run) + classical
    `C-mvalue-enetS`; 4096-col univariate prefilter; CLI split modes.
  - EWAS Atlas overlap screen on panel CpGs (validation only).
  - **Three-stage trait universe** on nine-pack full gene-linked (age/sex/tissue +
    disease/cancer): catalog → restrict → fold-safe reselect.
  - G1 comparison smoke vs 12b all-genes **not run**.
- **Results (CPU dry-run ATS, 2026-09-09):** see
  [`reports/inspection/stage0_7g_prime_matched_probe/analysis.md`](../../reports/inspection/stage0_7g_prime_matched_probe/analysis.md)
  + `cpu_panel_overview.json`.
  - **3 panels** (folds 0–2); traits **age / sex / tissue**; universe **65 536**
    CpGs → **2 658** genes in assignment.
  - Panel sizes ~**33–35k** CpGs → ~**1.2–1.3k** genes / fold; **717** genes
    shared across all folds (union **1 776**).
  - Classical mean: age r **0.894**, sex AUROC **0.891**, tissue bal-acc **0.437**.
  - Atlas: age/sex panel CpGs **enriched** vs curated associations; tissue Atlas
    in-universe too sparse for a meaningful enrichment call.
- **Results (trait universe):** fill from
  `reports/inspection/stage0_11_trait_universe/` when catalog + fold-reselect
  complete.
- **Verdict:** ATS CPU prep **usable** for G1 discussion; **not** a G1 lock.
  Trait-universe path is the scalable discovery design. Does **not** block
  N-light 65k OOF; neural Stage B still explicit-schedule.
## Policy (2026-09-09)

Milestone **10** locked the architecture finalists (**P2-G** cascade +
**N-light@64**). Milestone **11** does **not** block **N-light 65k** OOF.
It **does** participate in **GATE G1** gene-set choice before product cascade
(see [`milestone-12b-full-gene-panel.md`](milestone-12b-full-gene-panel.md) §2c).

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
  (**CPU dry-run 3/3 present**, `n_repeats=1`; production wants 5)
- Matched `C-mvalue-enetS` (**3/3 folds present**) / finalist-on-S (`N-cascade-S`,
  light) / optional fusion — neural **still open**
- Overview report:
  `reports/inspection/stage0_7g_prime_matched_probe/analysis.md`
- `direct_cpg.zarr` when direct loci exist (`n_direct > 0`)
- For **G1**: written gene-set verdict (adopt M11 panel vs all represented genes)

**Does not block:** Milestone **12 N-light** 65k validation.

**12b gene-set option:** a fold-selected gene list (this milestone / ADR 0012:
discovery CpGs → seed genes → all linked CpGs of those genes) is a valid
DeepRVAT-style training-set alternative to “all ~19.6k represented genes” in
[`milestone-12b-full-gene-panel.md`](milestone-12b-full-gene-panel.md) §2c.
Still deferred for implementation; in scope for GATE G1 before product cascade.

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

# Classical enetS after panels exist (partial --folds merges into existing JSON)
uv run python -u scripts/run_7g_prime_stage_b.py --device cpu --classical-only --folds 0,1,2

# Or: wait for fold-0 panel → panels 1–2 → classical all folds (BLAS threads capped)
PANEL_REPEATS=1 bash scripts/run_stage_b_cpu_chain.sh
```

Cap OpenBLAS when selecting on wide matrices:
`OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8` (set in the chain script).
Stability selection univariate-prefilters to **4096** columns before the enet
grid (same screen as gene-seed panels); meta records `n_cols_after_prefilter`.

`--panels-only` and `--classical-only` are mutually exclusive. Panels-only
writes `fold_panels/manifest.json` (merges fold rows on re-run) and does **not**
regenerate the full Stage B report. Classical-only loads existing
`fold_*_panel.json` and writes `per_arm/C-mvalue-enetS.json` (merges fold rows).
`--panel-repeats` sets stability-selection repeats (default 5; dry-run `1` is
recorded in the panel JSON/manifest).
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
