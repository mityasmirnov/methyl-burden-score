# Data population — what we have

Stage 0 open data under `$MBS_DATA_ROOT` / release `deepmat-data-v1`.

**Live numbers** (regenerate after catalog refresh):

```bash
source scripts/activate_data_environment.sh
uv run python scripts/write_data_population_inventory.py
# → reports/inspection/deepmat_data_v1/data_population_inventory.{md,json}
```

Narrative + HTTP trees: [`EWAS_DATA.md`](EWAS_DATA.md). Phenotype contracts:
[`EWAS_METADATA.md`](EWAS_METADATA.md). Catalog layout: [`DATA_CATALOG.md`](DATA_CATALOG.md).

## Snapshot (2026-09-28)

| Layer | Count |
|-------|------:|
| Catalog samples | **185 851** |
| Catalog studies | **1 860** |
| Catalogued assay files | **187 324** |
| GEO metadata GSM (parquet) | **170 338** |
| Census ∩ catalog overlap | **169 367** |
| Census not in catalog | **11 267** (Hub-empty / partial / ID aliases — not local download debt) |
| Empty EWAS_db dirs | **129** (remote also empty; prior refill lists fully applied) |
| Null `study.platform_id` | **42** studies / **~9.8k** samples (**39** mixed-platform — correctly null) |

### Platforms (studies)

| Platform | Studies |
|----------|--------:|
| HM450 | 1 162 |
| EPIC | 636 |
| EPICv2 | 20 |
| (null, mostly mixed) | 42 |

Live gap report: `make residual-gaps` →
`reports/inspection/deepmat_data_v1/residual_gaps.md`.
Assay probe-count platform fill: `make residual-gaps-apply-platforms`
(updates `catalog/tables/study.parquet`; DuckDB syncs on next
`catalog-refresh-release`).

### Phenotypes (distinct GSM in `sample_phenotype`)

| Trait | GSM | Notes |
|-------|----:|-------|
| tissue | 157 909 | mapped categorical |
| sex | 117 671 | |
| age | 100 260 | ~90 910 with numeric |
| disease | 90 877 | case 73 711 / control 10 441 / unknown 6 930 |
| bmi | 33 652 | |
| cancer | 15 247 | case 12 682 / control 688 / unknown 1 877 |
| blood | 3 402 | Hub pack |
| brain | 1 997 | Hub pack |
| ancestry | 1 380 | fairness / domain eval |

**Trait eligibility:** 47 rows; **24** core-eligible. Full gate table in the
inventory report. Sources mix Hub nine-pack sample-info, GEO brief-SOFT, EWAS
repository metadata, and Atlas study context.

### Top tissues (catalog)

| Tissue | GSM |
|--------|----:|
| whole blood | 41 934 |
| PBMC | 8 345 |
| leukocyte | 8 117 |
| brain | 8 045 |
| breast | 5 362 |
| peripheral blood | 3 308 |
| lung | 3 274 |
| cord blood | 3 205 |
| liver | 3 105 |
| placenta / kidney / saliva / … | see inventory |

### Assay mirror (`EWAS_db/`)

- Remote index: **1 989** studies (all dirs present on disk).
- Usable profiles: **1 860** dirs with `*.txt`, **187 324** sample files
  (**173 166** `GSM*.txt`).
- Empty dirs: **129** — Hub remote also has **0** `*.txt` (probed 2026-09-28);
  cannot refill from the mirror.
- Prior empty-dir refill lists (`ewas_db_refill_*.txt`) are **fully applied**.
- Commands: `make residual-gaps` / `make residual-gaps-apply-platforms` /
  `MBS_SKIP_ATLAS_SEED=1 make catalog-refresh-release` /
  `make data-population-inventory`.

### GEO labels (complete for catalog GSE)

Brief SOFT metadata (not betas): **170 338** GSM / **1 707** primary studies;
age ~65k, sex ~110k; series title on **1 718** GSE, `overall_design` on **271**.
See [`plans/improve-labels-study-context.md`](plans/improve-labels-study-context.md).
