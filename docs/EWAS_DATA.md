# EWAS Open Platform downloads

**Primary Stage 0 open methylation source.** Per
[ADR 0002](adr/0002-ewas-datahub-primary-source.md) and
[`STRATEGIC_PLAN.md`](STRATEGIC_PLAN.md), EWAS Data Hub supplies the default
pilot / open-scale matrices; EWAS Atlas supplies curated associations for later
validation. CpGCorpus remains optional (see [`CPGCORPUS_STAGE0.md`](CPGCORPUS_STAGE0.md)).

Stage 0 pulls **both** EWAS Atlas (association knowledge) and EWAS DataHub
(normalized methylation profiles). Raw files land under `$MBS_DATA_ROOT/raw/`.

Portal: https://ngdc.cncb.ac.cn/ewas/  
DataHub download page: https://ngdc.cncb.ac.cn/ewas/datahub/download  
Atlas downloads: https://ngdc.cncb.ac.cn/ewas/downloads

## Policy

- **Wanted: all public DataHub data** — both *All Data* (`EWAS_db`) and
  *Baseline Data* packs listed below.
- **Wanted: full Atlas batch exports** (associations, studies, cohorts,
  probe annotations, trait–trait relationships).
- Store under `data/raw/ewas_datahub/` and `data/raw/ewas_atlas/` only.
  Never under `vendor/`.
- Prefer resume-capable downloads (`wget -c`). Dedicated FTP clients
  (FileZilla) are recommended by NGDC for the large FTP trees.

## EWAS Atlas → `data/raw/ewas_atlas/`

Script: `scripts/download_ewas_atlas.sh` / `make download-ewas-atlas`

| File | Role |
|------|------|
| `EWAS_Atlas_associations.tsv` | curated EWAS associations |
| `EWAS_Atlas_studies.tsv` | studies |
| `EWAS_Atlas_cohorts.tsv` | cohorts |
| `EWAS_Atlas_probe_annotations.tsv` | probe annotations |
| `EWAS_trait_trait_logP.txt` | trait–trait relationships |

HTTP base used by the script: `https://download.cncb.ac.cn/ewas/`

## EWAS DataHub → `data/raw/ewas_datahub/`

Script: `scripts/download_ewas_datahub.sh` / `make download-ewas-datahub`

HTTP root (preferred on this host):
[https://download.cncb.ac.cn/ewas/datahub/](https://download.cncb.ac.cn/ewas/datahub/)

Remote index children:

```text
EWAS_db/        # All Data — per-study GSM*.txt beta files (~1989 studies)
add_ewas_db/    # supplemental (add_txt_450 / add_txt_850 / add_txt_935; may be empty)
download/       # Baseline packs (*_methylation_v1.zip, sample_*.zip, GMQN.zip)
```

Downloader: `scripts/download_ewas_datahub.sh` parses nginx HTML indexes and
`wget -c` each file (recursive wget alone is unreliable here because of
JS-enhanced listings / robots.txt). Use `--execute robots=off` is not enough
for `EWAS_db`; the script lists `href`s explicitly.

Inspection snapshot (hosts, layouts, local examples):
[`reports/inspection/raw_data_snapshot/summary.md`](../reports/inspection/raw_data_snapshot/summary.md).

**Full size / schema / top-row inventory (for organizing ingest):**
[`reports/inspection/raw_inventory/summary.md`](../reports/inspection/raw_inventory/summary.md).

**Atlas small tables + Hub sample-info structure contracts:**
[`docs/EWAS_METADATA.md`](EWAS_METADATA.md) and
[`reports/inspection/ewas_metadata_structure/`](../reports/inspection/ewas_metadata_structure/)
(`mbs inspect ewas-metadata`).

Local layout after mirror:

```text
data/raw/ewas_datahub/
  EWAS_db/          # All Data only — provenance lane ewas_datahub_db
  add_ewas_db/
  download/         # Baseline packs only — lane ewas_datahub_baseline
  SOURCE.txt
```

Do **not** leave baseline zips as flat `ewas_datahub/*.zip`; everything under
`download/`. Atlas stays in `raw/ewas_atlas/`; CpGCorpus Arrow in
`raw/cpgcorpus/` (Stage 0 GSEs) with aborted-sync leftovers in
`raw/cpgcorpus/_partial_fullsync/`. Catalog seeds these as separate
`provenance_lane` / `source_system` values (see `sql/002_provenance_lanes.sql`).

### All Data (`EWAS_db/`)

NGDC also advertises FTP (FileZilla): `ftp://download.big.ac.cn/ewas/datahub/EWAS_db/`  
Prefer HTTP: `https://download.cncb.ac.cn/ewas/datahub/EWAS_db/`

Per-study layout: `{STUDY}/GSM*.txt` with `probe_id<TAB>beta` (no header).
**All 18 Stage 0 labeling GSEs are present** here.

### Baseline Data (`download/`)

HTTP: `https://download.cncb.ac.cn/ewas/datahub/download/`  
FTP: `ftp://download.big.ac.cn/ewas/datahub/download/`

Matrix packs stay under `download/`. **Sample-info zips** are downloaded, then
unzipped to `reports/inspection/ewas_datahub_samples/` (zips may be deleted
after extract). Prefer the unpacked `.txt` for export; see
[`EWAS_METADATA.md`](EWAS_METADATA.md).

| Description | Archive | Approx. size |
|-------------|---------|--------------|
| DNA methylation profiles of 31 organism parts | `tissue_methylation_v1.zip` | 7.7 GB |
| Sample information (tissue) | `sample_tissue_methylation_v1.zip` → `…/sample_tissue.txt` | 62 KB |
| DNA methylation profiles of 25 brain parts | `brain_methylation_v1.zip` | 2.77 GB |
| Sample information (brain) | `sample_brain_methylation_v1.zip` → `…/sample_brain.txt` | 27 KB |
| DNA methylation profiles of 25 blood cell types | `blood_methylation_v1.zip` | 4.86 GB |
| Sample information (blood) | `sample_blood_methylation_v1.zip` → `…/sample_blood.txt` | 42 KB |
| Male/female profiles in 24 tissues | `sex_methylation_v1.zip` | 4.33 GB |
| Sample information (sex) | `sample_sex_methylation_v1.zip` → `…/sample_sex.txt` | 38 KB |
| DNA methylation changes with age | `age_methylation_v1.zip` | 11.73 GB |
| Sample information (age) | `sample_age_methylation_v1.zip` → `…/sample_age.txt` | 89 KB |
| Six ancestry categories | `ancestry_category_methylation_v1.zip` | 1.96 GB |
| Sample information (ancestry) | `sample_ancestry_category_methylation_v1.zip` → **`sample_race.txt`** | 21 KB |
| DNA methylation changes with BMI | `bmi_methylation_v1.zip` | 3.06 GB |
| Sample information (BMI) | `sample_bmi_methylation_v1.zip` → `…/sample_bmi.txt` | 28 KB |
| Profiles of 39 cancers | `cancer_methylation_v1.zip` | 16.07 GB |
| Sample information (cancer) | `sample_cancer_methylation_v1.zip` → `…/sample_cancer.txt` | 117 KB |
| Profiles of 28 diseases | `disease_methylation_v1.zip` | 20.11 GB |
| Sample information (disease) | `sample_disease_methylation_v1.zip` → `…/sample_disease.txt` | 154 KB |
| GMQN reference materials | `GMQN.zip` | ~40 MB |

All nine sample-info families above are unpacked and covered by
`mbs inspect ewas-metadata` (as of 2026-08-06). Canonical Parquet exports live
under `$MBS_DATA_ROOT/canonical/phenotypes/{family}_sample_info.parquet`.

### Local download status (refresh via inventory)

Authoritative on-disk snapshot:
[`reports/inspection/raw_inventory/summary.md`](../reports/inspection/raw_inventory/summary.md)
(+ narrative in [`DATA_CATALOG.md`](DATA_CATALOG.md)).

| Asset | Status (2026-08-24) |
|-------|---------------------|
| All nine Hub **profile** zips (age…disease) | **Complete** (Zip OK; disease finished 2026-08-11 after CNCB drop/416 failures) |
| All nine Hub **sample-info** zips + Parquet | **Complete** |
| Atlas batch TSVs | **Complete** (~0.26 GiB) |
| EPICv2 manifests | **Complete** |
| `EWAS_db/` All Data | **In progress** — see ingest counts below (~990 GiB betas under `ewas_datahub/`) |
| Host disk | `/data` ~2 T free at last check — packs are fine; watch space as `EWAS_db` grows |

### EWAS_db ingest counts (what the release manifest reports)

`mbs catalog refresh-release` does a **shallow directory listing** of
`$MBS_DATA_ROOT/raw/ewas_datahub/EWAS_db/` (no beta reads, no per-file content
hashes). Snapshot from the 2026-08-24 `deepmat-data-v1` refresh:

| Field | Value | Meaning |
|-------|-------|---------|
| `n_local_studies` | **883** | Study dirs under `EWAS_db/` that contain at least one `*.txt` GSM beta file |
| `n_local_gsm` | **87 153** | Total `*.txt` GSM files found across those studies |
| `advertised_n` | **1989** | Study count advertised by the remote EWAS_db index (or the last successful `--fetch-remote-index` fetch; default constant is 1989) |
| `mirror_complete` | **false** | `n_local_studies < advertised_n` — All-Data mirror is still downloading |

Empty study dirs (wget not started or not finished) are skipped and do **not**
count toward `n_local_studies`. More study dirs may exist on disk than 883; only
dirs with `GSM*.txt` are ingested. Hub nine-pack phenotypes are a separate SoT
and are already complete — EWAS_db completeness is **not** a Milestone 7A/7B gate
([ADR 0007](adr/0007-crossfit-prerequisites.md)).

After additional `EWAS_db/{GSE}/` study dirs appear (or in-progress downloads finish),
re-run the versioned catalog release — it rescans the tree and upserts studies
without rebuilding Hub phenotypes from scratch:

```bash
uv run mbs catalog refresh-release
# or: make catalog-refresh-release
uv run mbs catalog validate-release
uv run mbs catalog phenotype-census
```

Release layout: `$MBS_DATA_ROOT/canonical/releases/deepmat-data-v1/`.
Census reports: `reports/inspection/deepmat_data_v1/`.
See [`plans/milestone-7a-harmonized-release.md`](plans/milestone-7a-harmonized-release.md).

Disease resume helper (size + EOCD gate; survives refused / bogus 416):

```bash
nohup bash scripts/download_disease_pack_resilient.sh \
  > "$MBS_ARTIFACT_ROOT/logs/downloads/disease_resilient_nohup.log" 2>&1 &
```

EWAS_db resume:

```bash
nohup bash scripts/download_ewas_datahub.sh EWAS_db \
  > "$MBS_ARTIFACT_ROOT/logs/downloads/ewas_datahub_EWAS_db.log" 2>&1 &
```

Phenotype sample counts (unique GSM): age 8,374; tissue 5,323; sex 2,978;
blood 3,402; brain 1,997; BMI 2,070; ancestry 1,380; cancer 10,101;
disease 12,218. See catalog for study Ns and matrix conversions.
## Matrix convert (what it means)

**Convert** turns downloaded Hub **profile packs** (large ZIP of per-sample
betas) into a **canonical matrix store** the trainer can stream:

| Step | Input | Output |
|------|--------|--------|
| Select samples | `*_sample_info.parquet` + study filter / `--all-studies` / `max_per_study` | Sample list |
| Stream betas | Hub profile ZIP (or EWAS_db study folder) | `betas.zarr` (samples × loci) |
| Sidecars | Probe map + phenotypes | `locus_index.parquet`, `sample_index.parquet`, `sample_phenotypes.parquet`, `matrix_manifest.json` |

CLI:

```bash
# Single EWAS_db study (pilot)
uv run mbs matrix convert --study-id GSE35069 --platform-id HM450 --verify

# Hub phenotype pack → study-subset or full matrix
uv run mbs matrix convert-pack \
  --phenotype-family age \
  --study-ids GSE51032,GSE56105 \
  --matrix-id matrix-hub-age-studyholdout-v2 \
  --max-per-study 100 \
  --platform-id HM450
```

Convert is **I/O + rewrite**, not model training. Caps (`max_per_study`) exist
only to keep Stage 0 holdouts small; Milestone 5d drops caps (`--all-studies`).
Scripts: `scripts/convert_hub_pack_subsets.sh`, `scripts/convert_hub_full_packs.sh`.

### Supplemental (`add_ewas_db/`)

HTTP: `https://download.cncb.ac.cn/ewas/datahub/add_ewas_db/`  
Contains `add_txt_450/`, `add_txt_850/`, `add_txt_935/` additions to the All Data corpus.

## GMQN (DataHub normalization)

DataHub profiles are prepared with Gaussian Mixture Quantile Normalization
(GMQN), a reference-based method that removes technical variation at signal
intensity level for 450K and EPIC/850K arrays (type I Gaussian mixture
rescaling, then type II via BMIQ; reference study GSE105018).

Citation: *GMQN: A Reference-Based Method for Correcting Batch Effects and
Probe Bias in HumanMethylation BeadChip.* Front. Genet. 2022.
[PMID=35069703](https://pubmed.ncbi.nlm.nih.gov/35069703/)

## Commands

```bash
source scripts/activate_data_environment.sh

# Atlas batch exports
nohup bash scripts/download_ewas_atlas.sh \
  > "$MBS_ARTIFACT_ROOT/logs/downloads/ewas_atlas.log" 2>&1 &

# DataHub HTTP trees (EWAS_db + add_ewas_db + download)
nohup bash scripts/download_ewas_datahub.sh all \
  > "$MBS_ARTIFACT_ROOT/logs/downloads/ewas_datahub.log" 2>&1 &

# Or one tree at a time:
# bash scripts/download_ewas_datahub.sh EWAS_db
# bash scripts/download_ewas_datahub.sh download
# bash scripts/download_ewas_datahub.sh add_ewas_db

# Single phenotype family (profile + sample-info) — Milestone 5b wave-1 order:
# age → tissue → disease
make download-ewas-family FAMILY=age
make download-ewas-family FAMILY=tissue
make download-ewas-family FAMILY=disease
# Background example:
# nohup make download-ewas-family FAMILY=age \
#   > "$MBS_ARTIFACT_ROOT/logs/downloads/ewas_family_age.log" 2>&1 &

# Export sample-info → canonical Parquet (prefers unpacked .txt, else zip)
make export-ewas-sample-info FAMILY=tissue
# Structure profile for Atlas small tables + sample packs:
uv run mbs inspect ewas-metadata
# R fallback when only .RData is present:
# Rscript scripts/export_ewas_sample_info.R tissue \
#   "$MBS_DATA_ROOT/raw/ewas_datahub/download/sample_tissue_methylation_v1.zip" \
#   "$MBS_DATA_ROOT/canonical/phenotypes/tissue_sample_info.parquet"

# Phenotype registry (git): configs/data/phenotype_registry.yaml
# Checksums after download: $MBS_DATA_ROOT/canonical/registries/download_checksums.parquet

# Single EWAS_db study (preferred for Stage 0 pilot; ~60 GSM files for GSE35069)
make download-ewas-study STUDY=GSE35069
# or: bash scripts/download_ewas_datahub_study.sh GSE35069

# Convert Hub study → canonical matrix store
uv run mbs matrix convert \
  --study-id GSE35069 \
  --platform-id HM450 \
  --verify
```

Inspect with DuckDB / sanitized reports only; do not recursively index raw
archives in the coding agent.
