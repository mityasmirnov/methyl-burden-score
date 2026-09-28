# Plan: Age-first seed-mask screen (7G′ Stage B blocker)

> **Canonical milestone: 9c**. Historical 7G′ age-primary seed-mask. See [`MILESTONE_INDEX.md`](MILESTONE_INDEX.md).

Status (2026-09-06): **`done` — seed-masking not adopted.**

## Question / Approaches / Results / Verdict

| Field | Content |
|-------|---------|
| **Question** | Does age-primary phenotype masking to fold-safe seed genes help vs dense all-gene cascade on ATS? |
| **Approaches** | G0 all-gene control; G1 seed heads only; G2 expanded seed-gene CpGs + masks; G3 matched-random; C0/C2 classical enet comparators. Fold 0, seeds {42,43}, 40 epochs after bug fixes. |
| **Results** | Mean G0 age MAE **17.0**, tissue F1 **0.233**, sex AUROC **0.788**. Masked arms (G1–G3) age MAE **21–25**, tissue F1 **≤0.10**, sex AUROC mostly **0.51–0.72**. C0 classical age MAE **8.94** still leads. ADR 0012 panel: age discovery **1024** → expanded **9595** CpGs (seed frac ≈9%). |
| **Verdict** | **Do not adopt seed-gene masking** for pretrained MBS/RBS age-primary training. G0 beats every masked variant; G1/G2/G3 are not separable once training converges. Proceed with dense cascade / N-light paths (Milestones 10–12). |

Normative report: [`reports/inspection/stage0_7g_prime_seed_mask/analysis.md`](../../reports/inspection/stage0_7g_prime_seed_mask/analysis.md).
Panel audit: [`panel_audit.md`](../../reports/inspection/stage0_7g_prime_seed_mask/panel_audit.md) (`ok_for_seed_mask_gpu: true`).

### Bugs fixed before the conclusive run

1. Classical `KeyError: 'tissues'` — thread `class_names` into `_phenotype_arrays`.
2. Silent `learning_rate` threading bug (`lr=1e-2` instead of `1e-3`) + missing gradient clipping → near-total cascade collapse at 15 epochs; fixed, budget raised 15→40.

### Fold-0 panel (ADR 0012 empirical confirmation)

| Trait | Prefilter | Discovery CpGs | Seed genes | Unique expanded | Seed frac |
|-------|----------:|---------------:|-----------:|----------------:|----------:|
| age | 4096 | 1024 | 256 | 9595 | 0.089 |
| tissue | 4096 | 45 | 41 | 2778 | 0.016 |
| sex | 4096 | 44 | 50 | 2356 | 0.019 |
| sex_autosome | 4096 | 44 | 50 | 2356 | 0.019 |

Gene union across age/tissue/sex: **331**. Expanded CpG union: **11 785**.
Age used univariate top-k fallback (`sparsity_ok: false` for age only);
tissue/sex are sparse. `graph_content_hash` non-null; `panel_hash`
`ef6cd307513f28e3c45021f85c3a4d0c`.

### Final GPU metrics (fold 0, 40 epochs)

| Arm | Seed | Age MAE | Tissue F1 | Sex AUROC | Best epoch |
|-----|-----:|--------:|----------:|----------:|-----------:|
| G0 | 42 | 17.66 | 0.227 | 0.790 | 40/40 |
| G0 | 43 | 16.38 | 0.238 | 0.786 | 40/40 |
| G1 | 42 | 21.48 | 0.098 | 0.512 | 38/40 |
| G1 | 43 | 25.13 | 0.058 | 0.515 | 39/40 |
| G2 | 42 | 22.12 | 0.070 | 0.565 | 40/40 |
| G2 | 43 | 21.48 | 0.074 | 0.715 | 39/40 |
| G3 | 42 | 24.26 | 0.074 | 0.583 | 23/40 |
| G3 | 43 | 25.20 | 0.000 | 0.489 | 2/40 (collapsed) |
| C0 | — | 8.94 | 0.367 | 0.854 | — |
| C2 | — | 10.61 | 0.352 | 0.856 | — |

## Done already (scaffolding + audit)

| Item | Evidence |
|------|----------|
| ADR 0011 / 0012 + runner/YAML | `configs/experiment/stage0_7g_prime_seed_mask.yaml`, `scripts/run_7g_prime_seed_mask.py` |
| Fold-0 `internal_fold` panel | `seed_panels/fold_0/` — hashes above |
| Provenance audit green | `panel_audit.md` |
| GPU screen conclusive | `analysis.md` § GPU screen results |

## Follow-ons (out of this milestone)

- Platform-agnostic robustness / CpG dropout → **GATE G3** (Milestone 12c), not reopened here.
- Fold-selected panels → **Milestone 11**.
- Do **not** regenerate fold-0 panels unless graph hash or selection code changes.

Normative: [ADR 0011](../adr/0011-seed-gene-sources.md),
[ADR 0012](../adr/0012-seed-gene-discovery-vs-deployment-input.md),
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
| Encoder | P2-G max/max, `explicit_only` |
| Selection | `validation_age_mae` → tissue F1 → sex AUROC |
| Loss | `lambda_age=1.0`, `lambda_tissue=0.3`, `lambda_sex=0.1` |
| Grid | Fold 0, seeds {42, 43}, K=256 |
| Masks | Age / tissue / sex `SeedMaskedLinearHead`; fail if &lt;32 genes |
| Discovery CpGs | Fold-safe association ranks genes only (ADR 0012) |
| G2 / C2 CpGs | **All** explicit gene-linked CpGs of selected genes |
| Traits | Config-driven; ATS default age / tissue / sex |
| Adoption | **Rejected** — dense G0 preferred |

## Screening grid

| Arm | CpGs | Head genes | Question |
|-----|------|------------|----------|
| G0 | All gene-linked | All genes | Age-primary all-gene control |
| G1 | All gene-linked | Trait seed masks | Supervision masking? |
| G2 | All expanded CpGs of seed genes | Trait seed masks | Input filtering? |
| G3 | Matched random | Matched random masks | Biological specificity? |
| C0 | G0 CpGs | Ridge/enet | Classical all-gene |
| C2 | Exact G2 expanded CpGs | Ridge/enet | Fair classical seed comparator |

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
| BMI | **blocked** on ATS | Eligible on `matrix-hub-bmi-full-v1` |
| disease / cancer | **blocked** | unknown ≠ control |
| blood / brain subtraits | **blocked** | After ontology / label quality |

## Seed sources (ADR 0011)

| Source | Leakage |
|--------|---------|
| `external_clean` | Fixed prior — **not** fold-fitted |
| `internal_fold` | Fold-fitted (outer train only) |
| `hybrid_fold` | Fold-safe if combined inside train |

## Artifacts

- Config: `configs/experiment/stage0_7g_prime_seed_mask.yaml`
- Runner: `scripts/run_7g_prime_seed_mask.py`
- Report: `reports/inspection/stage0_7g_prime_seed_mask/`
- Panels: `seed_panel.json` + gene/locus parquet (hashed)

## Later: platform-agnostic robustness (GATE G3)

Deployment aggregates observed eligible CpGs per gene (coordinate + build;
never impute absent as zero). Covered under Milestone **12c / G3**, not 9c.

## Non-goals

- Re-adopting seed masks for age-primary ATS training after this negative result
- BMI / disease / cancer heads on ATS
- Treating age univariate fallback as sparse elastic-net seeds
