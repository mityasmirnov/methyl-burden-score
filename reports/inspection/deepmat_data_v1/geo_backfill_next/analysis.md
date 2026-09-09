# GEO backfill next (Wave 2b) — complete

**Date:** 2026-09-08

## Method

Switched fetch from FTP `*_family.soft.gz` (embeds full methylation tables,
often 1–2 GB+) to NCBI Accession Display **brief** SOFT
(`acc.cgi?targ=all&view=brief&form=text`) — metadata only, typically hundreds
of KB. Cached as `cache/geo/{GSE}/{GSE}_family.brief.soft.gz`. Default CLI:
`--view brief` (FTP `full` remains as fallback).

Stream-parse still skips any embedded tables if a full SOFT is used.

## Results

| Item | Value |
|------|------:|
| Studies requested | 100 |
| Download/parse ok | **100** |
| Failures | **0** |
| GEO parquet GSM (all fetched studies) | **111 195** (305 GSE) |
| Next-100 GSM | **19 363** |
| Next-100 age / sex / tissue mapped | **7 855** / **11 945** / **6 111** |
| Series overall_design after enrich | **267** / 1 718 (was 111) |

See `fetch_status.json` and plan `docs/plans/improve-labels-study-context.md`.
