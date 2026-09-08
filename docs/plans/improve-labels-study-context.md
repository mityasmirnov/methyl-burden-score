# Improve labels and study context beyond Hub packs

**Status:** Wave 1 + Wave 2a **done** (2026-09-07); Wave 2b list next  
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

## Wave 2b — next audited GEO crawl list

```bash
make write-geo-next-gse-list
# → configs/data/geo_backfill_next_gse.txt + reports/.../geo_backfill_next/
```

## Wave 2 (remaining, gated)

| Item | Gate |
|------|------|
| Fetch next GSE list (SOFT download) | Wave 2b list reviewed |
| Disease/cancer from tumor source_name + case/control pairing report | Per-study audit green |
| GEO-dev matrix convert + age/tissue train arms | Milestone 10/11 allow |
| Platform-null residual fill (`geo_platform_gap_priority_gse.txt`) | Ops parallel |
| overall_design for more studies | Needs family SOFT cache (download) |

## Non-goals

- Treating Atlas cohort fields as sample labels
- ComBat / study ID as encoder features
- Claiming GEO disease ready for training
- Full EWAS_db SOFT crawl in this wave

## Commands

```bash
source scripts/activate_data_environment.sh
# Rebuild phenotypes from cached SOFT (no FTP)
uv run python scripts/fetch_geo_sample_metadata.py \
  --studies-file configs/data/geo_backfill_batch50_gse.txt \
  --from-cache-only
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
