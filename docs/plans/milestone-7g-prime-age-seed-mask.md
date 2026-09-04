# Plan: Age-first seed-mask screen (7G′ Stage B blocker)

Status: **blocked on 16-ep GPU only** — fold-0 seed-panel auditability gate is
met (`graph_content_hash` set, stability diagnostics ≠ prefilter width for
tissue/sex, age uses documented `univariate_prefilter_top_k` when SGD coefs
explode, `sex_autosome` control with zero X/Y seeds, overlap + G3 quality in
`analysis.md` / `panel_audit.md`). Do **not** launch
`scripts/run_7g_prime_seed_mask.py --device cuda` while the 16-ep queue owns
GPU 0. `--reuse-panels` is allowed only for panels with non-null
`graph_content_hash` (runner enforces). Stage B CpG-panel GPU stays blocked
until this screen and typed-RBS diagnostics land.

## Done already (2026-09-04)

| Item | Evidence |
|------|----------|
| ADR 0011 / 0012 + runner/YAML scaffolding | `configs/experiment/stage0_7g_prime_seed_mask.yaml`, `scripts/run_7g_prime_seed_mask.py` |
| Fold-0 `internal_fold` panel regenerated | `seed_panels/fold_0/seed_panel.json` — `panel_hash` `ef6cd307…`, `graph_content_hash` `7ee70c55…` |
| Provenance / sparsity audit green | `panel_audit.md` → `ok_for_seed_mask_gpu: true` |
| Tissue/sex discovery sparse | 45 / 44 discovery CpGs (`sparsity_ok`); age fallback documented |
| `sex_autosome` control | Present; `n_sex_chrom_seed_cpgs=0` |
| Overlap + G3 matching quality | In `analysis.md` / panel JSON (gene union 331; G3 exact CpG-count match ≈90.6%) |
| Atlas catalog (non-blocking) | `sql/013_*`, `association_catalog.py` |

## Next steps

1. Wait for matched 16-ep promotion to finish and refresh
   `promotion_decision.json` (see [`milestone-7g-prime-16ep-promotion.md`](milestone-7g-prime-16ep-promotion.md)).
2. Clear `scratch/SEED_MASK_GPU_BLOCKED.txt` only when `next_gate` unlocks
   age-primary seed-mask (or equivalent).
3. Launch CUDA screen:  
   `uv run python -u scripts/run_7g_prime_seed_mask.py --device cuda --reuse-panels`  
   (fold 0, seeds {42,43}, arms G0/G1/G2/G3/C0/C2).
4. Write seed-mask `summary.json` / analysis; only then consider Stage B
   fold-panel GPU (`run_7g_prime_stage_b.py` still blocked).

Do **not** regenerate fold-0 panels unless the graph hash or selection code
changes; do **not** `--reuse-panels` on the stale null-hash copy under
`fold_0.stale-null-graph-hash/` / `fold_0.pre-topk-fallback/`.

Normative: [ADR 0011](../adr/0011-seed-gene-sources.md) (seed sources),
[ADR 0012](../adr/0012-seed-gene-discovery-vs-deployment-input.md)
(discovery CpGs vs deployment input),
[ADR 0010](../adr/0010-gene-allocation-policy.md).
Parents: [`milestone-7g-prime-matched-probe-lightweight.md`](milestone-7g-prime-matched-probe-lightweight.md),
[`milestone-7g-prime-pre-stage-b.md`](milestone-7g-prime-pre-stage-b.md).

## Goal

Freeze **`P2-G` as the current reference, not a final lock.** Run an
age-primary, phenotype-masked seed-gene screen on **`internal_fold`** panels
with the G0/G1/G2/G3/C0/C2 decomposition.

## Locked decisions

| Choice | Decision |
|--------|----------|
| First seed source | `internal_fold` only (DeepRVAT-faithful) |
| Atlas `external_clean` / `hybrid_fold` | Parallel catalog track; not a GPU gate |
| Encoder | P2-G max/max 15 ep, `explicit_only` |
| Selection | `validation_age_mae` → tissue F1 → sex AUROC |
| Loss | `lambda_age=1.0`, `lambda_tissue=0.3`, `lambda_sex=0.1` |
| Grid | Fold 0, seeds {42, 43}, K=256 then maybe 512 |
| Masks | Age / tissue / sex `SeedMaskedLinearHead`; fail if &lt;32 genes |
| Discovery CpGs | Fold-safe association ranks genes only (ADR 0012) |
| G2 / C2 CpGs | **All** explicit gene-linked CpGs of selected genes |
| Traits | Config-driven; ATS default age / tissue / sex |

## Screening grid

| Arm | CpGs | Head genes | Question |
|-----|------|------------|----------|
| G0 | All gene-linked | All genes | Age-primary all-gene control |
| G1 | All gene-linked | Trait seed masks | Supervision masking? |
| G2 | All expanded CpGs of seed genes | Trait seed masks | Input filtering? |
| G3 | Matched random | Matched random masks | Biological specificity? |
| C0 | G0 CpGs | Ridge/enet | Classical all-gene |
| C2 | Exact G2 expanded CpGs | Ridge/enet | Fair classical seed comparator |

Discovery CpGs (stability / prefilter survivors, often 4,096) are **not** the
G2 input panel. G2 and C2 train on sibling-enriched gene CpGs; `is_seed_cpg`
marks which loci came from discovery.

## Required panel report fields (per trait)

| Field | Meaning |
|-------|---------|
| `n_discovery_cpgs` | CpGs surviving stability selection |
| `n_seed_genes` | selected genes |
| `n_expanded_gene_cpg_edges` | all selected gene–CpG edges |
| `n_unique_expanded_gene_cpgs` | unique CpGs used by G2/C2 |
| `n_multigene_cpgs` | CpGs attached to multiple selected genes |
| `seed_fraction_of_expanded` | unique discovery / unique expanded |

## Trait catalog

| Trait | Status on ATS seed-mask | Notes |
|-------|-------------------------|-------|
| age | **active** (primary) | Core continuous |
| tissue | **active** (secondary) | Core multiclass |
| sex | **active** (auxiliary) + `sex_autosome` control | Autosomal sensitivity |
| BMI | **blocked** on ATS | Eligible on `matrix-hub-bmi-full-v1` (2,070 / 25); fails ≥1k bar on ATS age-pack BMI |
| disease / cancer | **blocked** | Need documented multi-study cases+controls; unknown ≠ control |
| blood / brain subtraits | **blocked** | After ontology / label quality |

Experiment YAML drives the active list (`seed_panel.traits`). Unknown ids or
traits without label arrays fail closed. Heads stay age/tissue/sex until a
later eligible-trait change.

## Seed sources (ADR 0011)

| Source | Leakage |
|--------|---------|
| `external_clean` | Fixed prior (Atlas − overlapping studies) — **not** fold-fitted |
| `internal_fold` | Fold-fitted (outer train only) |
| `hybrid_fold` | Fold-safe if combined inside train |

## Artifacts

- Config: `configs/experiment/stage0_7g_prime_seed_mask.yaml`
- Runner: `scripts/run_7g_prime_seed_mask.py`
- Report: `reports/inspection/stage0_7g_prime_seed_mask/`
- Panels: `seed_panel.json` + gene/locus parquet (hashed)

## Later: platform-agnostic robustness (not this screen)

Deployment aggregates whatever eligible CpGs a gene has on 450K, EPIC, or
ONT (coordinate + build; pool observed only; never impute absent as zero).
Reuse presence-aware paths in
`mbs.training.transparent_baselines.presence_aware_means` and observed-edge
assembly in `mbs.training.features`. Before calling the encoder
platform-agnostic: EPIC→450K / heterogeneous-coverage dropout, score
stability vs gene coverage, min-coverage / low-confidence flags, and
fold-fitted platform transforms (do not feed ONT frequency into an array
M-value model as identical measurements). **No training code in this plan.**

## Non-goals

- Launching Stage B CpG-panel GPU or Milestone 7
- BMI / smoking / disease / cancer heads (or joining BMI onto ATS)
- 450K↔EPIC dropout, ONT transforms, coverage-confidence tensors
- Per-tissue-class masks in the first grid
- Blocking on full Atlas DuckDB ingest
- Treating age univariate fallback as sparse elastic-net seeds (age
  `ranking_fallback=univariate_prefilter_top_k` when SGD coefs explode)
