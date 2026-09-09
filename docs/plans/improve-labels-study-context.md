# Improve labels and study context beyond Hub packs

**Status:** Wave 1–2c **done** (catalog-wide GEO sample metadata, 2026-09-09)  
**Parent:** [`data-infrastructure-improvements.md`](data-infrastructure-improvements.md) §2  
**Related:** [`geo-metadata-backfill-ewas-db.md`](geo-metadata-backfill-ewas-db.md),
[`geo-enriched-training-release.md`](geo-enriched-training-release.md),
[`EWAS_METADATA.md`](../EWAS_METADATA.md)

## Scope and acceptance

**Problem:** Beyond Hub nine-pack sample-info, most catalog GSMs are EWAS_db-only.
GEO batch-50 + Atlas enrichment exist but under-deliver:

| Gap | Evidence (2026-09-07) |
|-----|------------------------|
| Tissue labels missing | **26 242** / 56 201 GEO rows `tissue_map_status=empty` despite `source_name` often holding tissue (e.g. GSE197678 PBMC) — parquet built **before** source_name fallback |
| Tissue unmapped | **15 296** rows; many strings already aliasable (`PBL`, `Blood DNA`, `DLPFC`) or are tumor→organ |
| Atlas study context | Only **182** / **1 763** catalog studies matched (~10%) |
| GEO study coverage | **111** studies in parquet vs **~1 700** EWAS_db studies |
| Disease/cancer | Global eligibility can look green; per-study case/control often unbalanced |

**Done when (this wave):**

1. Rebuild `geo_sample_metadata.parquet` from SOFT cache with source_name tissue
   fallback + expanded aliases; census GEO tissue mapped count rises materially.
2. Re-seed Atlas GSE↔ES map against **all catalog GSE** (not only GEO-parquet
   studies); report matched N in `study_atlas_enrichment.md`.
3. Catalog refresh writes updated census / eligibility / Atlas enrichment.
4. This plan updated with before/after numbers; training still **not** wired to
   GEO disease/cancer.

## Locked decisions

| Choice | Decision | Why |
|--------|----------|-----|
| Mutate ATS / Hub packs? | **No** | Frozen SoT |
| Tumor `source_name` → tissue | Map to **organ** tissue when clear (`brain tumor`→`brain`, `breast tumor`→`breast`); do **not** invent disease controls | Tissue ≠ diagnosis |
| Tumor names as disease labels | Separate follow-on (case-only without controls stay non-core) | DATA_CONTRACT |
| Encoder features from GEO text | **No** | Phenotypes / study metadata only |
| Larger GEO crawl | After this rebuild wave | Fix what we already fetched first |

## Wave 1 (this change set)

1. Expand `configs/data/geo_tissue_aliases.yaml` for high-count empty/unmapped strings
   that map safely to Hub ontology organs/cell types.
2. `uv run python scripts/fetch_geo_sample_metadata.py --studies-file … --from-cache-only`
   on batch-50 (+ any cached studies).
3. `make seed-atlas-gse-map` over catalog GSE list (existing script; ensure input
   is catalog studies, not only parquet).
4. `MBS_SKIP_ATLAS_SEED=1 make catalog-refresh-release` (or full refresh after seed).
5. Record Δ: mapped tissue N, Atlas matched N, top residual gaps.

## Wave 2a — GEO series study context (**done** 2026-09-07)

| Metric | Value |
|--------|------:|
| Catalog GSE with `metadata_json.geo` | **1 718** / 1 718 |
| With title / summary | **1 718** / **1 718** |
| With overall_design (SOFT header) | **111** (cached family SOFT only) |
| NCBI esummary | 1 607; esummary+SOFT 111 |

Artifacts: `canonical/phenotypes/geo_series_metadata.parquet`,
`reports/inspection/deepmat_data_v1/geo_series_enrichment.{json,md}`.
Merged on every `mbs catalog refresh-release` via `merge_geo_series_into_studies`.

```bash
make enrich-geo-series-metadata
MBS_SKIP_ATLAS_SEED=1 make catalog-refresh-release
```

## Wave 2b — next audited GEO crawl (**done** 2026-09-08)

```bash
make write-geo-next-gse-list
# → configs/data/geo_backfill_next_gse.txt (100 GSE)
uv run python scripts/fetch_geo_sample_metadata.py \
  --studies-file configs/data/geo_backfill_next_gse.txt --view brief
make enrich-geo-series-metadata
MBS_SKIP_ATLAS_SEED=1 make catalog-refresh-release
```

**List (2026-09-08):** **440** candidates (≥50 GSM); **100** selected.

**Fetch (same day):** **100/100** ok, **0** failures via NCBI Accession Display
``view=brief`` (metadata-only SOFT; KB–MB vs multi-GB FTP ``*_family.soft.gz``).
Default fetch mode is now ``--view brief`` (FTP ``full`` kept as fallback).

| Metric | Value |
|--------|------:|
| GEO parquet GSM | **111 195** (305 studies) |
| Next-100 GSM added | **19 363** |
| Next-100 with age / sex / tissue mapped | **7 855** / **11 945** / **6 111** |
| Species (parquet) | human **110 838**; non-human **357** (quarantined at merge) |
| Series `overall_design` (SOFT brief) | **267** / 1 718 (was 111) |

Artifacts: `reports/inspection/deepmat_data_v1/geo_backfill_next/fetch_status.json`,
`cache/geo/*/GSE*_family.brief.soft.gz`.

## Wave 2c — remaining catalog GSE (**done** 2026-09-09)

Prior remain-all crawl (2026-09-08) brought parquet to ~1 639 studies / 162k GSM.
Final gap: **79** catalog GSE missing as primary `study_id`:

```bash
# configs/data/geo_backfill_catalog_remain_gse.txt
uv run python scripts/fetch_geo_sample_metadata.py \
  --studies-file configs/data/geo_backfill_catalog_remain_gse.txt --view brief
make enrich-geo-series-metadata
MBS_SKIP_ATLAS_SEED=1 make catalog-refresh-release
```

| Metric | Value |
|--------|------:|
| Fetch | **79/79** ok, **0** failures |
| GEO parquet | **170 338** GSM / **1 707** primary studies |
| Catalog GSE | **1 718** |
| Primary-study gap | **11** related/sub-series (100% GSM overlap via `study_ids`) |
| age / sex / tissue mapped | **64 727** / **110 430** / **70 244** |
| Series `overall_design` | **271** / 1 718 |

GEO **sample metadata** for catalog GSE is complete. Remaining data-population
work is EWAS_db **assay** mirror (`mirror_complete=false`, ~1 695/1 989 studies).

Report: `reports/inspection/deepmat_data_v1/geo_backfill_remain/analysis.md`.

## Wave 2 (remaining, gated)

| Item | Gate |
|------|------|
| Disease/cancer from tumor source_name + case/control pairing report | Per-study audit green |
| GEO-dev matrix convert + age/tissue train arms | Milestone 10/11 allow |
| Platform-null residual fill (`geo_platform_gap_priority_gse.txt`) | Ops parallel |
| EWAS_db assay mirror completion / failure retry | Download ops |

## Non-goals

- Treating Atlas cohort fields as sample labels
- ComBat / study ID as encoder features
- Claiming GEO disease ready for training
- Treating GEO metadata complete as EWAS_db assay mirror complete

## Commands

```bash
source scripts/activate_data_environment.sh
# Preferred: metadata-only brief SOFT (no methylation tables)
uv run python scripts/fetch_geo_sample_metadata.py \
  --studies-file configs/data/geo_backfill_next_gse.txt --view brief
# Rebuild phenotypes from cache only
uv run python scripts/fetch_geo_sample_metadata.py \
  --studies-file configs/data/geo_backfill_batch50_gse.txt \
  --from-cache-only
make enrich-geo-series-metadata
make seed-atlas-gse-map
make catalog-refresh-release
```

## Acceptance evidence

- `reports/inspection/deepmat_data_v1/geo_backfill_batch/` (rebuild summary)
- `reports/inspection/deepmat_data_v1/study_atlas_enrichment.md` (matched N ↑)
- `reports/inspection/deepmat_data_v1/census.md` (GEO phenotype rows)
- Before/after table in this file (fill after run)

## Before / after (Wave 1 — 2026-09-07)

| Metric | Before | After |
|--------|-------:|------:|
| GEO parquet `tissue` mapped | 14 663 | **43 886** |
| GEO parquet empty tissue | 26 242 | **0** |
| Catalog GEO tissue rows (mapped only) | ~21k (mixed raw) | **36 861** (0 `nan`) |
| Catalog GEO age / sex GSM | — | **21 940** / **32 083** |
| Top new tissue classes | — | whole blood 12 950; `brain - tumor` 10 252; leukocyte 4 336; PBMC 3 331 |
| Atlas matched studies | 182 / 1 763 | **182 / 1 763** (PMID overlap ceiling; 1 490 GSE have PubMed, 182 in Atlas) |

**Ops notes:** Full SOFT re-parse of all 111 cached family files blew past ~48 GB RAM — use
`uv run python scripts/remap_geo_tissue.py [--force]` instead. Catalog writes
**mapped** tissue only (unmapped raw stays on parquet). Atlas match rate is limited
by GEO↔Atlas PMID overlap, not seed coverage; study context beyond Atlas needs
GEO series metadata (Wave 2).

```bash
uv run python scripts/remap_geo_tissue.py --force
MBS_SKIP_ATLAS_SEED=1 make catalog-refresh-release
```
