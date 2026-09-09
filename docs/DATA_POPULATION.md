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

## Snapshot (2026-09-09)

| Layer | Count |
|-------|------:|
| Catalog samples | **173 076** |
| Catalog studies | **1 763** |
| Phenotype rows | **673 786** |
| Catalogued assay files | **170 641** |
| GEO metadata GSM (parquet) | **170 338** |
| Trait eligibility rows | **47** (24 core-eligible) |

### Platforms (samples via study.platform_id)

| Platform | Samples |
|----------|--------:|
| HM450 | 105 639 |
| EPIC | 55 901 |
| (null) | 10 415 |
| EPICv2 | 1 121 |

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
- Usable profiles: **~1 682** dirs with `GSM*.txt`, **~159.5k** GSM files,
  **~1 728** dirs with any `*.txt` (rising; ~261 empty mid-refill).
- Gap: recoverable plan = **162** GSE (~15.9k GSM; ~129 still empty mid-pass) +
  **3** remaining empty non-GSM (HCMI-CMDC, TARGET-ALL-P3, TARGET-AML) —
  `ewas_db_empty_refill_plan.json`. Parallel refill:
  `EWAS_REFILL_JOBS` × `EWAS_REFILL_FILE_JOBS`.
- `mirror_complete` stays false until recoverable empties are filled.
- Commands: `make refill-ewas-db-empty` /
  `LIST=configs/data/ewas_db_refill_nongsm.txt make refill-ewas-db-empty` →
  then `MBS_SKIP_ATLAS_SEED=1 make catalog-refresh-release` and
  `make data-population-inventory`.

### GEO labels (complete for catalog GSE)

Brief SOFT metadata (not betas): **170 338** GSM / **1 707** primary studies;
age ~65k, sex ~110k; series title on **1 718** GSE, `overall_design` on **271**.
See [`plans/improve-labels-study-context.md`](plans/improve-labels-study-context.md).
