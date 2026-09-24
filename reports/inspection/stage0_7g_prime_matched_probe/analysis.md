# Milestone 11 CPU panel overview (ATS 3-fold dry-run)

**As of:** 2026-09-09  
**Machine-readable:** [`cpu_panel_overview.json`](cpu_panel_overview.json)  
**Panels:** [`fold_panels/`](fold_panels/) · **Classical:** [`per_arm/C-mvalue-enetS.json`](per_arm/C-mvalue-enetS.json)  
**Orphan census:** [`orphan_rbs_census.md`](orphan_rbs_census.md)

## What was done

CPU-only Milestone **11** prep on ATS (`hub-ats-7e-3fold-v1`), **65 536-prefix**
HM450 columns, graph-v2:

1. Honest Stage A lock refresh (provisional **P2-G**, next gate = scale review).
2. Stage B CLI: `--panels-only` / `--classical-only` / `--folds` / `--panel-repeats`.
3. Univariate **4096-col** prefilter before stability SGD (was multi-day on 65k).
4. **3 fold panels** written (`n_repeats=1` dry-run; production should use **5**).
5. Classical **`C-mvalue-enetS`** on each fold’s panel (all 3 folds).
6. Orphan RBS census + 10c trait-hygiene policy (no phenotype rewrite).
7. EWAS Atlas overlap screen (validation only; not used for selection).

**Not done:** neural Stage B (`N-cascade-S`, N-light-on-S, fusion), production
`n_repeats=5` panels, GATE **G1** gene-set verdict vs all ~19.6k genes.

## Scope

| Item | Value |
|------|------:|
| Traits | age, sex, tissue |
| Folds / panels tested | **3** (one panel per outer fold) |
| Matrix | `matrix-hub-age-tissue-sex-full-v1` |
| Split | `hub-ats-7e-3fold-v1` |
| CpG universe scored | **65 536** prefix |
| Genes in assignment (65k) | **2 658** |
| Selector | study-grouped multitask enet stability |
| Prefilter | top **4 096** univariate cols / trait |
| Stability repeats | **1** (dry-run) |

## Panel sizes (CpGs → genes)

Graph expansion: seed CpGs → all gene-linked CpGs of those genes (ADR 0012 style).

| Fold | Train N | Seed CpGs | Panel CpGs | Seed genes | Panel genes | Seeds by trait (age/sex/tissue) |
|-----:|--------:|----------:|-----------:|-----------:|------------:|--------------------------------|
| 0 | 8 163 | 3 410 | 32 731 | 1 035 | 1 207 | 3 333 / 43 / 34 |
| 1 | 7 207 | 3 408 | 34 784 | 1 106 | 1 272 | 3 333 / 34 / 50 |
| 2 | 7 631 | 3 419 | 35 182 | 1 150 | 1 302 | 3 333 / 45 / 45 |

**Cross-fold:**

| Quantity | Intersection (all 3) | Union | Mean pairwise Jaccard |
|----------|---------------------:|-------:|------------------------:|
| Panel CpGs | 23 097 | 44 373 | ~0.66 |
| Panel genes | **717** | **1 776** | ~0.57 |
| Seed CpGs | — | — | ~0.25 (0–1 weak; 1–2 higher) |

Age dominates seed quota (`max_seeds/3 ≈ 3333`); frequency gate failed for age
(`sparsity_ok: false`, `n_passing_min_frequency: 0`) so age seeds are the
prefilter/top-k fallback head. Sex/tissue sparsity gates passed (~34–50 cols).

## Classical `C-mvalue-enetS` (held-out fold)

| Fold | Age Pearson r | Sex AUROC | Tissue bal-acc | Tissue macro-F1 |
|-----:|--------------:|----------:|---------------:|----------------:|
| 0 | 0.894 | 0.852 | 0.422 | 0.368 |
| 1 | 0.877 | 0.928 | 0.466 | 0.395 |
| 2 | 0.912 | 0.892 | 0.424 | 0.393 |
| **mean** | **0.894** | **0.891** | **0.437** | **0.385** |

Age/sex are strong on the fold-selected sparse panel; tissue multi-class remains
hard (expected for ~40 scored classes on ATS).

## EWAS Atlas correlation (validation)

Source: `data/raw/ewas_atlas/EWAS_Atlas_associations.tsv` (latin-1). Trait labels
matched by substring/exact rules (see JSON). Universe = 65k `probe_id`s.
Hypergeometric enrichment of **panel CpGs** vs Atlas probes in that universe:

| Trait | Atlas probes in 65k | Mean panel overlap | Mean fraction of panel | Mean hypergeom p |
|-------|--------------------:|-------------------:|-----------------------:|-----------------:|
| age | 3 197 | ~2 248 | **6.6%** | **≪ 10⁻⁷⁵** (enriched) |
| sex | 1 377 | ~792 | **2.3%** | **~4×10⁻⁴** (enriched) |
| tissue | 124 | ~64 | **0.19%** | ~0.55 (not enriched) |

**Read:** age and sex fold-panels recover Atlas-curated loci far above chance
inside the 65k prefix. Tissue Atlas coverage in-universe is tiny (124 probes),
so the non-significant tissue enrichment is **not** evidence against the panel —
Atlas tissue labels are sparse relative to Hub CE tissue.

## Orphan RBS (related product note)

At the same 65k × graph-v2 × `explicit_only`: **1 953** orphan regions
(**861** singleton / **1 092** multi-CpG). Multi-CpG orphans are the fusion-
qualified set; see orphan census.

## Caveats

- Dry-run **`n_repeats=1`** — production panels need **5** before G1 adoption.
- Age seed list is quota-saturated / frequency-gate fallback — interpret seed
  counts as “prefilter head,” not stable multitask frequency hits.
- Atlas rules are heuristic string matches; not a genome-wide EWAS re-analysis.
- Neural Stage B arms not run; classical ≠ cascade product score.

## Pointers

- Plan: [`docs/plans/milestone-11-fold-selected-panel.md`](../../../docs/plans/milestone-11-fold-selected-panel.md)
- TODO §11: [`docs/TODO_PIPELINE.md`](../../../docs/TODO_PIPELINE.md)
- Runner: `scripts/run_7g_prime_stage_b.py` / `scripts/run_stage_b_cpu_chain.sh`
