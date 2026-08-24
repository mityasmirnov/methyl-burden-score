# Raw data inventory (refreshed)

Inspected: **2026-08-24T10:32:27.865942+00:00**. Machine-readable: [`summary.json`](summary.json).

Shallow sizes for Hub profile packs, Atlas, EPICv2 manifests, and Stage 0 CpGCorpus GSEs. Does not recurse into every `EWAS_db` GSM file.

## Totals under `data/raw/`

| Tree | Bytes (approx.) |
|------|----------------:|
| `ewas_datahub/` | 990.03 GiB |
| `cpgcorpus/` | 6.72 GiB |
| `manifests/` | 445.1 MiB |
| `ewas_atlas/` | 268.7 MiB |
| `raw_total/` | 997.51 GiB |

EWAS_db study directories: **1049**

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
| disease | `disease_methylation_v1.zip` | 20.11 | 20.11 GiB | True | `ok` |

GMQN.zip: 40.1 MiB status=`ok`

## Hub sample-info zips (phenotypes, not betas)

| Family | File | On-disk |
|--------|------|--------:|
| age | `sample_age_methylation_v1.zip` | 227 KiB |
| tissue | `sample_tissue_methylation_v1.zip` | 191 KiB |
| sex | `sample_sex_methylation_v1.zip` | 112 KiB |
| blood | `sample_blood_methylation_v1.zip` | 91 KiB |
| brain | `sample_brain_methylation_v1.zip` | 48 KiB |
| bmi | `sample_bmi_methylation_v1.zip` | 103 KiB |
| ancestry | `sample_ancestry_category_methylation_v1.zip` | 44 KiB |
| cancer | `sample_cancer_methylation_v1.zip` | 485 KiB |
| disease | `sample_disease_methylation_v1.zip` | 394 KiB |

## Hub sample-info Parquet (`canonical/phenotypes/*_sample_info.parquet`)

Row counts can exceed unique GSM (duplicate rows in Hub R tables). Use **unique `sample_id`** as training N.

| Family | Rows | Unique GSM | Unique studies |
|--------|-----:|-----------:|---------------:|
| age | 8,374 | 8,374 | 143 |
| tissue | 5,323 | 5,323 | 258 |
| sex | 2,978 | 2,978 | 161 |
| blood | 3,402 | 3,402 | 161 |
| brain | 1,997 | 1,997 | 40 |
| bmi | 2,070 | 2,070 | 25 |
| ancestry | 1,380 | 1,380 | 21 |
| cancer | 10,841 | 10,101 | 225 |
| disease | 14,501 | 12,218 | 209 |

## EWAS_db All-Data tree (in progress)

Local study directories: **1049** / advertised remote **1989** (52.7% of study folders if the remote count is stable).
Per-GSM text files under `raw/ewas_datahub/EWAS_db/{GSE}/`. Resume: `bash scripts/download_ewas_datahub.sh EWAS_db`.

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
| GSE116992 | GPL13534 | 113.6 MiB |
| GSE116992 | GPL21145 | 288.7 MiB |
| GSE125367 | GPL21145 | 415.0 MiB |
| GSE35069 | GPL13534 | 312.4 MiB |

## Regenerate

```bash
uv run python scripts/write_raw_inventory_refresh.py
uv sync --extra analysis  # once, for matplotlib
uv run python scripts/write_pipeline_doc_figures.py
```
