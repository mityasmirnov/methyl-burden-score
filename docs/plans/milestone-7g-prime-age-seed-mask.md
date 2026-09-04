# Plan: Age-first seed-mask screen (7G′ Stage B blocker)

Status: **in progress** — next scientifically useful GPU run after Stage A
evidence. Stage B CpG-panel GPU (`C-mvalue-enetS` / `N-cascade-S`) stays
blocked until this screen and typed-RBS diagnostics land.

Normative: [ADR 0011](../adr/0011-seed-gene-sources.md),
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

## Screening grid

| Arm | CpGs | Head genes | Question |
|-----|------|------------|----------|
| G0 | All gene-linked | All genes | Age-primary all-gene control |
| G1 | All gene-linked | Trait seed masks | Supervision masking? |
| G2 | Union of seed-gene CpGs | Trait seed masks | Input filtering? |
| G3 | Matched random | Matched random masks | Biological specificity? |
| C0 | G0 CpGs | Ridge/enet | Classical all-gene |
| C2 | Exact G2 CpGs | Ridge/enet | Fair classical seed comparator |

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

## Non-goals

- Launching Stage B CpG-panel GPU or Milestone 7
- BMI / smoking / disease heads
- Per-tissue-class masks in the first grid
- Blocking on full Atlas DuckDB ingest
