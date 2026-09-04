# EWAS metadata contracts (Atlas small tables + Hub sample-info)

Short durable contracts for Cursor-visible EWAS tables. Full profiles live in
[`reports/inspection/ewas_metadata_structure/`](../reports/inspection/ewas_metadata_structure/)
(regenerate with `mbs inspect ewas-metadata`).

## Scope

| Lane | Files | Role |
|------|-------|------|
| `ewas_atlas` | `EWAS_Atlas_studies.tsv`, `EWAS_Atlas_cohorts.tsv`, `EWAS_trait_trait_logP.txt` | Association knowledge metadata / trait similarity |
| `ewas_datahub_baseline` | `reports/inspection/ewas_datahub_samples/sample_*/*.txt` | Sample-level phenotypes for Hub packs |

Out of scope here: large Atlas associations / probe annotations, matrix zips,
`EWAS_db` beta text.

## Parse recipes

| Source | How to read |
|--------|-------------|
| Atlas studies / cohorts | Tab-separated TSV with header (`pd.read_csv(..., sep="\t", encoding="latin-1")`; rare extra-tab rows skipped) |
| Atlas trait×trait | TSV; column 0 = trait name; remaining columns = square logP matrix |
| DataHub sample-info `.txt` | R `write.table` (space-separated, double quotes, row names). Use `mbs.registry.sample_info.read_r_style_table` (`index_col=0` dropped) |

## Join keys

| Table | Keys |
|-------|------|
| Atlas studies | `study_ID` (ES…), `PMID` |
| Atlas cohorts | `cohort_ID`, `study_ID` → studies |
| Hub sample-info | `sample_id` (GSM…), `project_id` (usually GSE…) |

**Do not** join Atlas `study_ID` to Hub `project_id` by raw string equality —
different namespaces (ES* vs GSE*). Prefer PMID / curated maps when linking
Atlas knowledge to Hub profiles.

Nine families currently unpacked and profiled: age, ancestry (`sample_race.txt`),
blood, bmi, brain, cancer, disease, sex, tissue.

## Family → primary phenotype column

Used by `mbs.registry.sample_info.FAMILY_VALUE_COLUMN` when exporting Parquet:

| Family | Column |
|--------|--------|
| age | `age` |
| tissue | `tissue` |
| disease | `disease` |
| cancer | `disease` |
| blood | `cell_component` |
| brain | `tissue` |
| sex | `sex` |
| ancestry | `race` (file: `sample_race.txt` inside `sample_ancestry_category_methylation_v1/`) |
| bmi | `bmi` |

## Export paths

Sample-info export prefers unpacked extracts, then download zips:

```text
reports/inspection/ewas_datahub_samples/sample_{family}_methylation_v1/sample_{family}.txt
# ancestry exception: …/sample_ancestry_category_methylation_v1/sample_race.txt
# else
$MBS_DATA_ROOT/raw/ewas_datahub/download/sample_{family}_methylation_v1.zip
→ $MBS_DATA_ROOT/canonical/phenotypes/{family}_sample_info.parquet
```

## Study-level Atlas enrichment (catalog)

On `mbs catalog refresh-release`, the release builder attaches **study-level**
Atlas context for external stratification (not sample phenotype labels):

- Parquet: `deepmat-data-v1/catalog/tables/study_atlas_enrichment.parquet`
- DuckDB: `study_atlas_enrichment`, view `v_study_atlas_enrichment`
- Report: `reports/inspection/deepmat_data_v1/study_atlas_enrichment.{json,md}`
- Also merged into `study.metadata_json` under `atlas_enrichment`

Join rules (never raw `GSE*` = `ES*` equality):

1. Optional curated map [`configs/data/atlas_gse_es_map.tsv`](../configs/data/atlas_gse_es_map.tsv)
   (`gse_id`, `atlas_study_id`, `pmid`, `source`)
2. PMID bridge when the map supplies `pmid` without `atlas_study_id`
3. Direct when catalog `study_id` is already an Atlas `ES*` accession

Fields include cohort count, total Atlas sample size, tissues, cohort
descriptions (disease area), platforms, ancestries, and Atlas trait names.

```bash
source scripts/activate_data_environment.sh
make catalog-refresh-release
```

Populate the curated map as GSE↔publication links are verified (GEO, papers).
Runs automatically before every `make catalog-refresh-release` (including the
EWAS_db post-download hook). Manual re-seed:

```bash
source scripts/activate_data_environment.sh
make seed-atlas-gse-map          # or: make catalog-refresh-release
# skip NCBI: MBS_SKIP_ATLAS_SEED=1 make catalog-refresh-release
```

Sample-level GEO backfill for EWAS_db-only GSM is a **separate catalog lane**
from Atlas enrichment. Atlas stays study-level. GEO writes GSM phenotypes with
`source_family=geo_metadata_backfill`. Hub sample-info still wins on overlapping
GSM. Operator brief, join keys, GPL map, and audit:
[`plans/geo-metadata-backfill-ewas-db.md`](plans/geo-metadata-backfill-ewas-db.md).
Pre-scale / batch-50:
[`plans/geo-metadata-backfill-pre-scale.md`](plans/geo-metadata-backfill-pre-scale.md).
Training release design (not built; do not mutate ATS):
[`plans/geo-enriched-training-release.md`](plans/geo-enriched-training-release.md).

```bash
source scripts/activate_data_environment.sh
make fetch-geo-sample-metadata          # pilot 15 (re-validate with current code)
# Batch expand only after repaired-pilot audit is accepted:
# make fetch-geo-sample-metadata-batch
MBS_SKIP_ATLAS_SEED=1 make catalog-refresh-release
# Clean Δ recipe: MBS_SKIP_GEO_BACKFILL=1 refresh → assert 0 GEO rows → merge
```

Reports: `reports/inspection/deepmat_data_v1/geo_backfill_pilot/`
(`summary.*`, `validation.*`, `fetch_status.json`).

Age is unit-aware (years/months/weeks/days); tissue maps through Hub ontology
aliases; diagnosis-only disease text is **not** a case/control label.

```bash
source scripts/activate_data_environment.sh
uv run mbs inspect ewas-metadata
make export-ewas-sample-info FAMILY=age
```

## Related

- Downloads: [`EWAS_DATA.md`](EWAS_DATA.md)
- Inspection guide: [`DATA_INSPECTION.md`](DATA_INSPECTION.md)
- Registry: `configs/data/phenotype_registry.yaml`
- Plan: [`plans/ewas-metadata-structure.md`](plans/ewas-metadata-structure.md)
- Pipeline gate: [`TODO_PIPELINE.md`](TODO_PIPELINE.md) §5b′
