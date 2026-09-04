# GEO metadata backfill — fixes before scaling

**Status:** fixes 1–4 implemented in code (2026-09-04); scale still blocked until
verified on pilot rebuild + remaining checklist items (disease cases, full
per-GSE status, training release)  
**Parent:** [`geo-metadata-backfill-ewas-db.md`](geo-metadata-backfill-ewas-db.md)
(pilot **done** / audited) · [`data-infrastructure-improvements.md`](data-infrastructure-improvements.md) §2  
**Related:** [`DATA_CONTRACT.md`](../DATA_CONTRACT.md), [`EWAS_METADATA.md`](../EWAS_METADATA.md),
[`TODO_PIPELINE.md`](../TODO_PIPELINE.md) (current gate = 7G′ promotion screen)

## Why this exists

The 15-GSE pilot proved fetch → parquet → Hub-wins merge → census. It is **not**
ready for a full EWAS_db crawl or for training on GEO phenotypes. The gaps below
would silently corrupt age, tissue classes, multi-GSE GSM identity, disease
eligibility, and audit deltas if we scaled first.

## Recommended order (binding)

1. **Finish remaining matched 16-epoch promotion arms** (7G′ current gate).
2. **Correct GEO backfill code** (this plan § Fixes 1–4 minimum; 5–6 as ready).
3. **Expand GEO** from 15 studies to a larger audited batch (status table per GSE).
4. **Eligibility report** by trait and study (`geo_metadata_backfill` + Hub).
5. **New GEO-enriched development cohort** (separate release id; not frozen ATS).
6. **Only then** add eligible GEO samples to age/tissue training and seed-gene
   selection — with Hub-only / GEO-only / Hub+GEO / metadata-only arms.

Do **not** silently enlarge ATS (`matrix-hub-age-tissue-sex-full-v1`) or retrain
frozen v0.1.

## Done when (pre-scale gate)

| # | Fix | Acceptance |
|---|-----|------------|
| 1 | Age unit-aware | Years/months/weeks/days parsed explicitly; original string retained; unit tests for `"6 months"`, `"120 days"`, bare years |
| 2 | Tissue ontology | GEO tissue → existing coarse ontology; unmapped/ambiguous reported; equivalent labels (whole blood / peripheral blood / blood) collapse when ontology says so |
| 3 | No silent GSM dedup | Multi-GSE GSM: metadata agreement check; sample–study membership persisted; conflicts counted, not `keep="first"` |
| 4 | Clean before/after audit | Incremental test starts from catalog with **zero** `geo_metadata_backfill` rows and reports exact phenotype delta |
| 5 | Disease extraction | Structured case/control study-by-study; still never invent controls from diagnosis-free samples; report case vs control counts |
| 6 | Per-GSE status table | Every crawled GSE has download status, SOFT checksum, GEO/EWAS_db GSM counts, phenotype counts, ambiguous fields, ontology coverage, conflict counts |
| 7 | Separate training release | New phenotype/matrix release id after QC; comparison arms Hub-only / GEO-only / Hub+GEO / metadata-only confounding |

Fixes **1–4** are the hard coding gate before any batch larger than the pilot.
**5–7** may land with the first expanded batch / cohort build.

## Fix 1 — Age unit-aware

**Today:** first numeric token → years (`"6 months"` → 6).

**Required:**

- Parse unit tokens: `year(s)|yr|y`, `month(s)|mo`, `week(s)|wk`, `day(s)|d`.
- Convert to years for `sample_phenotype.numeric_value` (age head); keep 0–120 year
  range after conversion; out-of-range / unknown unit → unknown / omit observed.
- Persist original string under `sample.metadata_json.geo.age_raw` (and optionally
  `age_unit`, `age_value_native`).
- Unit tests: bare `45`, `45 years`, `6 months` → 0.5, `120 days` → ~0.329,
  sentinel ages still rejected.

## Fix 2 — Harmonize tissue labels

**Today:** raw `tissue` strings; no ontology pass.

**Required:**

- Map GEO tissue through the existing tissue ontology path used for Hub
  (`mbs.training.phenotype_table` / `tissue_ontology*.yaml`), or a GEO-specific
  alias table that feeds the same class ids.
- Write `sample.tissue_ontology_id` when mapped; keep raw in
  `sample.tissue_raw` / `metadata_json.geo`.
- Report: mapped / unmapped / ambiguous counts; top unmapped strings.
- Training must not treat `whole blood`, `peripheral blood`, and `blood` as
  distinct classes when ontology collapses them.

## Fix 3 — Do not silently deduplicate GSMs

**Today:** `drop_duplicates(subset=["sample_id"], keep="first")` in
`scripts/fetch_geo_sample_metadata.py`.

**Required:**

- One GSM may appear in related GSE / superseries SOFT files.
- Compare parsed phenotype + key metadata fields across appearances.
- If agree → single phenotype row; persist **sample–study membership**
  (all GSE links), e.g. `sample.metadata_json.geo.study_ids` and/or a small
  membership table / parquet.
- If disagree → do not silently pick first; emit conflict record (counts +
  examples in per-GSE status / pilot report); leave conflicting traits unknown
  or omit until audited.
- Unit tests: same GSM, two GSE, agreeing vs conflicting sex/age.

## Fix 4 — Improve before/after audit

**Today:** re-refresh snapshots `census.json` that already contains GEO rows, so
the report shows ~32 967 GEO rows before and after.

**Required:**

- Explicit incremental test path:
  1. Refresh with `MBS_SKIP_GEO_BACKFILL=1` (or empty/missing parquet) → assert
     `count(*) WHERE source_family='geo_metadata_backfill' = 0`.
  2. Snapshot compact census.
  3. Enable parquet + merge → snapshot after.
  4. Report exact Δ phenotype rows / observed GSM / per family.
- Do not embed full census blobs; keep compact before/after + merge_stats.
- Document that unique-GSM movement can still be EWAS_db mirror growth.

## Fix 5 — Expand disease extraction carefully

**Today:** pilot disease rows are **controls only** (correct under conservative
rules; not a disease head).

**Required:**

- Study-by-study structured case extraction (explicit case/patient/disease-status
  language, GEO `sample_type` when present, audited mappings).
- **Never** infer that diagnosis-free or unlabeled samples are controls.
- Keep `unknown ≠ control` for eligibility.
- Report case/control/unknown per GSE; GEO disease remains non-core until cutoffs
  pass.

## Fix 6 — Scale by study with explicit status

Before any full crawl, persist **per-GSE** (parquet or JSONL under
`reports/inspection/.../geo_backfill/`):

| Field | Purpose |
|-------|---------|
| `study_id` | GSE |
| `download_status` / `error` | ok / fail / skip |
| `soft_sha256` | cache integrity |
| `n_geo_gsm` | samples in SOFT |
| `n_match_ewas_db` | intersection with local GSM txt |
| `n_hub_skipped` | Hub-wins |
| phenotype counts by id | age/sex/tissue/disease/cancer |
| ambiguous / unmapped fields | age unit fails, tissue unmapped, … |
| ontology mapping coverage | tissue mapped fraction |
| metadata conflict counts | multi-GSE GSM disagreements |

Makefile / CLI should be able to resume by study and skip already-ok GSE when
checksum matches.

## Fix 7 — Separate training release

After QC on an expanded audited batch:

- New release id (e.g. `deepmat-data-geo-dev-v1`), **not** overwriting
  `deepmat-data-v1` census semantics for frozen Hub packs and **not** mutating
  ATS `matrix-hub-age-tissue-sex-full-v1`.
- Derive phenotype tables / matrix membership for GEO-enriched samples.
- Comparison arms (documentation + eval configs):

  | Arm | Role |
  |-----|------|
  | Hub-only | Baseline |
  | GEO-only | GEO phenotype coverage / leakage |
  | Hub + GEO | Combined eligible labels |
  | Metadata-only confounding | Study/platform leakage alarm (7E′ style) |

## Non-goals (still)

- Wiring GEO into Milestone **7** OOF before 7G′ Stage A/B lock
- Treating Atlas cohort fields as sample labels
- ComBat / batch as encoder features
- Full CNCB EWAS_db crawl as a training gate

## Implementation checklist

- [x] Age unit parser + tests + `age_raw` in metadata_json
- [x] Tissue ontology pass + unmapped report
- [x] Multi-GSE GSM membership + conflict handling (remove silent `keep="first"`)
- [x] Zero-GEO baseline → merge incremental audit path
  (authoritative Δ = `n_phenotype_rows_added` / in-memory before_merge;
  dirty disk census flagged; operator recipe in report notes)
- [ ] Per-GSE status artifact schema + writer
  (partial: `fetch_status.json` with download/checksum/tissue_map/conflicts)
- [ ] Disease case extraction playbook (study batches; no invented controls)
- [ ] Design ADR or plan addendum for GEO-enriched training release id + arms
- [ ] Expand beyond 15 GSE only after fixes 1–4 green

## Code touchpoints (when implementing)

| Area | Likely files |
|------|----------------|
| Parse / eligibility | [`src/mbs/geo_metadata.py`](../../src/mbs/geo_metadata.py) |
| Fetch / dedup | [`scripts/fetch_geo_sample_metadata.py`](../../scripts/fetch_geo_sample_metadata.py) |
| Merge / census snapshot | [`src/mbs/release.py`](../../src/mbs/release.py) |
| Tissue ontology | [`src/mbs/training/phenotype_table.py`](../../src/mbs/training/phenotype_table.py) |
| Tests | [`tests/unit/test_geo_metadata.py`](../../tests/unit/test_geo_metadata.py) |
| Pilot / batch reports | `reports/inspection/deepmat_data_v1/geo_backfill*/` |
