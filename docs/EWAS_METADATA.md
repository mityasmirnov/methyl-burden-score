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
