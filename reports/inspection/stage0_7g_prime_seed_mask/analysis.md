# 7G′ age-primary seed-mask screen

Selection: validation **age MAE** primary; tissue macro-F1 secondary; sex AUROC
tertiary. Seed source: **`internal_fold`** ([ADR 0011](../../../docs/adr/0011-seed-gene-sources.md)).
P2-G topology is the **current reference, not a final lock**. Stage B CpG-panel
GPU stays blocked until this screen finishes.

Config: `configs/experiment/stage0_7g_prime_seed_mask.yaml` · runner:
`scripts/run_7g_prime_seed_mask.py`.

## Grid

| Arm | Question |
|-----|----------|
| G0 | Age-primary all-gene control (not tissue-heavy P2-G) |
| G1 | All gene-linked CpGs + trait seed head masks |
| G2 | Seed-gene CpGs + same masks |
| G3 | Matched-random genes/CpGs + matched masks |
| C0 | Classical enet on G0 CpGs |
| C2 | Classical enet on exact G2 CpGs |

Folds: `[0]`; seeds: `{42, 43}`; K=`256`; loss λ = age `1.0` / tissue `0.3` / sex `0.1`.

## Fold-0 `internal_fold` panels (landed)

Artifacts: `seed_panels/fold_0/seed_panel.{json,gene.parquet,locus.parquet}`.

| Trait | n genes | n seed CpGs (prefilter→enet) | n_runs |
|-------|--------:|-----------------------------:|-------:|
| age | 256 | 4096 | 36 |
| tissue | 256 | 4096 | 36 |
| sex | 256 | 4096 | 36 |

Constructor: univariate prefilter (top 4096) → study-grouped enet stability
(2×2 folds × α/l1 grid) → explicit-edge gene enrichment. Sex panel records
autosome-only flag for reporting.

## Metrics status

See `summary.json` (populated as arms finish).

**GPU 0 blocked** until the 7G′ 16-epoch promotion screen unlocks
(`scratch/SEED_MASK_GPU_BLOCKED.txt`; unlock =
`reports/inspection/stage0_7g_gene_only_probe/promotion_decision.json`).
Do **not** launch `scripts/run_7g_prime_seed_mask.py` on CUDA until then.
Fold-0 `internal_fold` panels are already on disk; restart with
`--reuse-panels --device cuda` after unlock.
