# Raw data inventory (refreshed)

Inspected: **2026-08-11T14:38:46.237445+00:00**. Machine-readable: [`summary.json`](summary.json).

Shallow sizes for Hub profile packs, Atlas, EPICv2 manifests, and Stage 0 CpGCorpus GSEs. Does not recurse into every `EWAS_db` GSM file.

## Totals under `data/raw/`

| Tree | Bytes (approx.) |
|------|----------------:|
| `ewas_datahub/` | 428.45 GiB |
| `cpgcorpus/` | 6.72 GiB |
| `manifests/` | 0.43 GiB |
| `ewas_atlas/` | 0.26 GiB |
| `raw_total/` | 435.93 GiB |

EWAS_db study directories: **457**

## Figures

![Raw tree sizes](figures/raw_tree_sizes.png)

![Hub pack sizes](figures/hub_pack_sizes.png)

![Hub sample counts](figures/hub_sample_counts.png)

## Hub profile packs (`raw/ewas_datahub/download/`)

| Family | File | Advertised GB | On-disk | Zip OK | Status |
|--------|------|--------------:|--------:|:------:|--------|
| age | `age_methylation_v1.zip` | 11.73 | 11.73 GiB | True | `ok` |
| tissue | `tissue_methylation_v1.zip` | 7.7 | 7.70 GiB | True | `ok` |
| sex | `sex_methylation_v1.zip` | 4.33 | 4.33 GiB | True | `ok` |
| blood | `blood_methylation_v1.zip` | 4.86 | 4.86 GiB | True | `ok` |
| brain | `brain_methylation_v1.zip` | 2.77 | 2.77 GiB | True | `ok` |
| bmi | `bmi_methylation_v1.zip` | 3.06 | 3.06 GiB | True | `ok` |
| ancestry | `ancestry_category_methylation_v1.zip` | 1.96 | 1.96 GiB | True | `ok` |
| cancer | `cancer_methylation_v1.zip` | 16.07 | 16.07 GiB | True | `ok` |
| disease | `disease_methylation_v1.zip` | 20.11 | 4.27 GiB | False | `bad_zip` |

GMQN.zip: 0.04 GiB status=`ok`

## EWAS Atlas files

| File | Bytes |
|------|------:|
| `EWAS_Atlas_associations.tsv` | 106,146,695 |
| `EWAS_Atlas_cohorts.tsv` | 397,573 |
| `EWAS_Atlas_probe_annotations.tsv` | 174,062,150 |
| `EWAS_Atlas_studies.tsv` | 153,681 |
| `EWAS_trait_trait_logP.txt` | 986,305 |

## EPICv2 manifests

| File | Bytes |
|------|------:|
| `EPICv2_reannotated_manifest_v3.0.csv.gz` | 466,768,830 |
| `SOURCE.txt` | 215 |

## CpGCorpus Stage 0 GSEs

| GSE | GPL | On-disk |
|-----|-----|--------:|
| GSE116992 | GPL13534 | 0.11 GiB |
| GSE116992 | GPL21145 | 0.28 GiB |
| GSE125367 | GPL21145 | 0.41 GiB |
| GSE35069 | GPL13534 | 0.31 GiB |

## Regenerate

```bash
uv run python scripts/write_raw_inventory_refresh.py
uv sync --extra analysis  # once, for matplotlib
uv run python scripts/write_pipeline_doc_figures.py
```
