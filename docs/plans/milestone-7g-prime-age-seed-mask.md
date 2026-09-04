# Plan: Age-first seed-mask screen (7G′ Stage B blocker)

Status: **blocked** — wait for (1) matched 16-epoch promotion decision rules
([`milestone-7g-prime-16ep-promotion.md`](milestone-7g-prime-16ep-promotion.md))
and (2) seed-panel provenance/sparsity audit (non-null `graph_content_hash`,
stability-selection diagnostics distinct from the 4096 univariate prefilter,
autosome-only sex control, trait overlap report). Do **not** launch
`scripts/run_7g_prime_seed_mask.py --device cuda` while the 16-ep queue owns
GPU 0, and do **not** `--reuse-panels` until the fold-0 panel is regenerated
with those fixes. Stage B CpG-panel GPU stays blocked until this screen and
typed-RBS diagnostics land.

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
- Relieving the sparsity audit (`4096 == prefilter` still fails GPU go/no-go)
