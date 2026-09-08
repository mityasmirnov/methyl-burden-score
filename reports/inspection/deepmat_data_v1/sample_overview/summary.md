# Sample overview (Hub + EWAS_db)

Generated: `2026-09-08T11:06:40Z`

Rows: **173076** (one per catalog `sample_id`).

## ID convention

- `sample_id` = **GSM** (sample)
- `study_id` / `gse_id` = **GSE** (series/study)
- `geo_gpl_id` = **GPL** (GEO platform)

## Lane membership

- Hub baseline: 34234
- EWAS_db: 169176
- Hub and EWAS_db: 30334
- GEO SOFT fetched: 90751

## Platform (after inference)

| Platform | N |
|----------|--:|
| HM450 | 109903 |
| EPIC | 56711 |
| EPICv2 | 6415 |
| unknown | 47 |

## Species (all rows filled)

| Status | N |
|--------|--:|
| human | 172721 |
| non_human | 355 |

## Hub pack membership

| Pack | N |
|------|--:|
| age | 8374 |
| tissue | 5323 |
| sex | 2978 |
| disease | 12218 |
| cancer | 10101 |
| blood | 3402 |
| brain | 1997 |
| bmi | 2070 |
| ancestry | 1380 |

## Phenotype coverage (non-null)

| Field | N |
|-------|--:|
| age_years | 58054 |
| sex | 82833 |
| tissue | 167756 |
| bmi | 2070 |
| ancestry_label | 1380 |
| brain_label | 1997 |
| cell_component | 3280 |
| pubmed_ids | 159483 |

## Disease pack (nine-pack Hub)

- Pack membership: **12218**
- Cases (patients): **5264** (43.1% of pack; 15.4% of nine-pack)
- Controls: **6930**
- Adjacent normal: **24**
- Distinct diagnosis labels (cases): **28**

| Disease | Cases | % of cases |
|---------|------:|----------:|
| Alzheimer's disease | 945 | 17.95 |
| schizophrenia | 536 | 10.18 |
| systemic lupus erythematosus | 341 | 6.48 |
| Parkinson's disease | 333 | 6.33 |
| ulcerative colitis | 258 | 4.9 |
| multiple sclerosis | 228 | 4.33 |
| rheumatoid arthritis | 225 | 4.27 |
| psoriasis | 211 | 4.01 |
| stroke | 204 | 3.88 |
| childhood asthma | 202 | 3.84 |
| intellectual disability and congenital anomalies | 200 | 3.8 |
| Crohn's disease | 197 | 3.74 |
| asthma | 194 | 3.69 |
| preeclampsia | 179 | 3.4 |
| Huntington's disease | 170 | 3.23 |
| systemic insulin resistance | 115 | 2.18 |
| respiratory allergy | 104 | 1.98 |
| type 2 diabetes | 92 | 1.75 |
| panic disorder | 89 | 1.69 |
| Silver Russell syndrome | 79 | 1.5 |
| Graves' disease | 73 | 1.39 |
| spina bifida | 59 | 1.12 |
| Down syndrome | 54 | 1.03 |
| autism spectrum disorder | 48 | 0.91 |
| Sjogren's syndrome | 48 | 0.91 |
| Kabuki syndrome | 40 | 0.76 |
| nephrogenic rest | 22 | 0.42 |
| systemic sclerosis | 18 | 0.34 |

Disease names come from Hub disease-pack phenotype_value on case rows; controls have sample_type=control and usually no diagnosis string. Train only on label_status in {case, control}.


## Notes

- Hub phenotype labels win; GEO fills blanks.
- Unknown platforms inferred from DataHub metadata, GPL maps, `_935k` hints, and EWAS_db probe counts.
- Species filled for every row (GEO SOFT > series taxon > Hub/project priors > EWAS assumed human); provenance in `species_source`.
- `cpgpt_sample_embedding_status=unavailable` (locus embeddings only).
- Artifacts: `data/canonical/phenotypes/sample_overview_hub_geo_v1.{parquet,csv.gz}`
- PDF dashboard: `sample_overview_dashboard.pdf`
- Canvas (project copy): `sample-overview-hub-geo.canvas.tsx`
