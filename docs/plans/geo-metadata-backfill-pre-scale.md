# GEO metadata backfill — fixes before scaling

**Status:** repaired **15-GSE pilot re-validated** (clean zero→merge Δ, 2026-09-04).
Batch expansion **blocked** until this audit is accepted. Training release
**designed only** ([`geo-enriched-training-release.md`](geo-enriched-training-release.md)).  
**Parent:** [`geo-metadata-backfill-ewas-db.md`](geo-metadata-backfill-ewas-db.md)
· [`data-infrastructure-improvements.md`](data-infrastructure-improvements.md) §2  
**Related:** [`DATA_CONTRACT.md`](../DATA_CONTRACT.md), [`EWAS_METADATA.md`](../EWAS_METADATA.md),
[`TODO_PIPELINE.md`](../TODO_PIPELINE.md)

## Why this exists

The original committed pilot summary was generated **before** unit-aware age,
tissue ontology, GSM conflict handling, and clean Δ audit — it proved the old
pipeline, not the repaired one. GEO still improves the **catalog only** (training
heads read Hub pack Parquet). Disease/cancer from GEO are **not** training-ready
in this pilot (disease controls-only; cancer cases without matched GEO controls /
eligibility). **Do not** silently enlarge frozen ATS.

## Recommended order (binding)

1. **Finish remaining matched 16-epoch promotion arms** (7G′ current gate; parallel OK).
2. **Correct GEO backfill code** (fixes 1–4) — **done** in tree.
3. **Re-validate the 15-study pilot with the new code** ← **done** (2026-09-04)
   - Rebuilt pilot parquet from cache.
   - Clean before/after: zero `geo_metadata_backfill` → merge → exact Δ (**33961**).
   - Tissue mapped/unmapped **per study** (`validation.md`).
   - Age distributions after unit conversion (`age_unit` all `years` in this slice).
   - Multi-GSE conflicts: **0** in pilot.
   - Per-GSE `fetch_status.json` + `validation.{json,md}`.
4. **Only then** expand to a larger audited GSE batch (e.g. batch-50) — **pending acceptance**.
5. Eligibility report by trait and study.
6. **Separate immutable** GEO-enriched training release
   (`deepmat-data-geo-dev-v1`) — samples with methylation **and** acceptable
   labels; Hub-only / GEO-only / Hub+GEO / metadata-only arms. Not ATS.
7. Only then add eligible GEO samples to age/tissue training / seed-gene work.

## Done when (pre-scale gate)

| # | Fix | Acceptance |
|---|-----|------------|
| 1 | Age unit-aware | Years/months/weeks/days parsed; `age_raw` retained; unit tests |
| 2 | Tissue ontology | Aliases + Hub map; unmapped/ambiguous reported |
| 3 | No silent GSM dedup | Membership + conflicts; not `keep="first"` |
| 4 | Clean before/after audit | Start from **zero** GEO rows; exact phenotype Δ |
| 4b | **Pilot re-validation report** | New `geo_backfill_pilot` summary after fixes; per-study tissue/age/conflict/status |
| 5 | Disease extraction | Explicit case/control only; never invent controls; **not for training while cases=0** |
| 6 | Per-GSE status table | Download, checksum, counts, ontology, conflicts |
| 7 | Separate training release | New immutable release id; not frozen ATS |

Code for 1–4 and 5 is in tree. **4b is the gate before step 4 (batch expand).**

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
  (**re-validated 2026-09-04:** `MBS_SKIP_GEO_BACKFILL=1` → 0 GEO rows → merge →
  disk_before=0, dirty=false, authoritative Δ=33961)
- [x] Per-GSE status + repaired-pilot validation
  (`fetch_status.json`, `validation.{json,md}`: tissue/age by study, conflicts)
- [x] Disease case extraction (explicit only; diagnosis-alone omitted)
  **Pilot after repair:** disease=2058 controls; cancer=994 cases (sample-type);
  **neither** trait is core-eligible — **do not train disease/cancer on GEO yet**
- [x] Training release **design only**
  → [`geo-enriched-training-release.md`](geo-enriched-training-release.md)
- [ ] Expand beyond 15 GSE — **blocked until this repaired-pilot audit is accepted**
  (batch-50 list/cache may exist; do **not** merge until then)

## Code touchpoints (when implementing)

| Area | Likely files |
|------|----------------|
| Parse / eligibility | [`src/mbs/geo_metadata.py`](../../src/mbs/geo_metadata.py) |
| Fetch / dedup | [`scripts/fetch_geo_sample_metadata.py`](../../scripts/fetch_geo_sample_metadata.py) |
| Merge / census snapshot | [`src/mbs/release.py`](../../src/mbs/release.py) |
| Tissue ontology | [`src/mbs/training/phenotype_table.py`](../../src/mbs/training/phenotype_table.py) |
| Tests | [`tests/unit/test_geo_metadata.py`](../../tests/unit/test_geo_metadata.py) |
| Pilot / batch reports | `reports/inspection/deepmat_data_v1/geo_backfill*/` |
