# EWAS Open Platform downloads

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

### All Data

NGDC recommends a dedicated FTP client:

```text
ftp://download.big.ac.cn/ewas/datahub/EWAS_db/
```

Local mirror target:

```text
$MBS_DATA_ROOT/raw/ewas_datahub/EWAS_db/
```

The downloader mirrors this tree with `wget --continue --recursive` when FTP
is reachable. If FTP is blocked or stalls, use FileZilla against the same host
and path, writing into that directory.

### Baseline Data

NGDC FTP host (alternate to HTTP):

```text
ftp://download.big.ac.cn/ewas/datahub/download/
```

HTTP mirrors used by the script (`*_v1.zip` names on `download.cncb.ac.cn`):

| Description | Archive | Approx. size |
|-------------|---------|--------------|
| DNA methylation profiles of 31 organism parts | `tissue_methylation_v1.zip` | 7.7 GB |
| Sample information (tissue) | `sample_tissue_methylation_v1.zip` | 62 KB |
| DNA methylation profiles of 25 brain parts | `brain_methylation_v1.zip` | 2.77 GB |
| Sample information (brain) | `sample_brain_methylation_v1.zip` | 27 KB |
| DNA methylation profiles of 25 blood cell types | `blood_methylation_v1.zip` | 4.86 GB |
| Sample information (blood) | `sample_blood_methylation_v1.zip` | 42 KB |
| Male/female profiles in 24 tissues | `sex_methylation_v1.zip` | 4.33 GB |
| Sample information (sex) | `sample_sex_methylation_v1.zip` | 38 KB |
| DNA methylation changes with age | `age_methylation_v1.zip` | 11.73 GB |
| Sample information (age) | `sample_age_methylation_v1.zip` | 89 KB |
| Six ancestry categories | `ancestry_category_methylation_v1.zip` | 1.96 GB |
| Sample information (ancestry) | `sample_ancestry_category_methylation_v1.zip` | 21 KB |
| DNA methylation changes with BMI | `bmi_methylation_v1.zip` | 3.06 GB |
| Sample information (BMI) | `sample_bmi_methylation_v1.zip` | 28 KB |
| Profiles of 39 cancers | `cancer_methylation_v1.zip` | 16.07 GB |
| Sample information (cancer) | `sample_cancer_methylation_v1.zip` | 117 KB |
| Profiles of 28 diseases | `disease_methylation_v1.zip` | 20.11 GB |
| Sample information (disease) | `sample_disease_methylation_v1.zip` | 154 KB |

Local layout:

```text
data/raw/ewas_datahub/
  EWAS_db/                 # All Data FTP mirror
  *_methylation_v1.zip     # Baseline packs
  sample_*_v1.zip
  SOURCE.txt
```

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

# DataHub: baseline packs + All Data FTP mirror
nohup bash scripts/download_ewas_datahub.sh \
  > "$MBS_ARTIFACT_ROOT/logs/downloads/ewas_datahub.log" 2>&1 &
```

Inspect with DuckDB / sanitized reports only; do not recursively index raw
archives in the coding agent.
