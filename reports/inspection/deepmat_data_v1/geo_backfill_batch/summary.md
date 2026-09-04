# GEO backfill pilot summary

- Generated: `2026-09-04T15:17:01Z`
- GEO parquet GSM in: **45762**
- Catalog samples touched (EWAS_db-only): **38020**
- Hub-skipped GSM: **7077**
- GEO GSM not in catalog: **665**
- Phenotype rows added: **74971**
- Samples with ≥1 observed GEO phenotype: **36014**
- Samples with `metadata_json.geo`: **38020**
- Multi-study GSM (membership): **0**

## Invariants

- GEO phenotype rows on Hub pack-membership GSM: **0** (must be 0)
- Atlas blobs on `sample.metadata_json`: **0** (must be 0)
- `sample_source_membership` includes `geo_metadata_backfill`: **False** (must be false)

## Tissue ontology

- mapped / unmapped / ambiguous / empty: **11670** / **9451** / **0** / **16899**
- Unmapped examples: `sarcoma`, `Medulloblastoma`, `brain tumor`, `cerebellum`, `frontal cortex`, `PBL`, `peripheral blood lymphocytes`, `baby_venous_blood`, `mother_venous_blood`, `breast tumor`, `Blood DNA`, `White blood cell`, `LM`, `PR`, `meningioma`, `ependymal tumor`, `Umbilical cord blood buffy coat`

## Phenotypes by id

| phenotype_id | rows | observed | unique GSM |
| --- | ---: | ---: | ---: |
| `age` | 19233 | 19233 | 19233 |
| `cancer` | 1470 | 1470 | 1470 |
| `disease` | 5216 | 5216 | 5216 |
| `sex` | 27931 | 27931 | 27931 |
| `tissue` | 21121 | 21121 | 21121 |

## Label status

| phenotype_id | label_status | n |
| --- | --- | ---: |
| `age` | `observed` | 19233 |
| `cancer` | `case` | 1239 |
| `cancer` | `control` | 231 |
| `disease` | `case` | 1456 |
| `disease` | `control` | 3760 |
| `sex` | `observed` | 27931 |
| `tissue` | `observed` | 21121 |

## Eligibility (`source_family=geo_metadata_backfill`)

| phenotype_id | n | cases | controls | core | aux |
| --- | ---: | ---: | ---: | --- | --- |
| `age` | 19233 | None | None | True | True |
| `cancer` | 1470 | 1239 | 231 | True | True |
| `disease` | 5216 | 1456 | 3760 | True | True |
| `sex` | 27931 | None | None | False | True |
  - `sex` not core: sex is auxiliary biological / QC, not a core burden target
| `tissue` | 21121 | None | None | True | True |

## Per study (merge)

| study_id | samples touched | phenotype rows |
| --- | ---: | ---: |
| `GSE104210` | 675 | 1245 |
| `GSE105018` | 1230 | 1230 |
| `GSE109379` | 1104 | 1104 |
| `GSE111629` | 316 | 1060 |
| `GSE116339` | 679 | 2037 |
| `GSE117859` | 608 | 1824 |
| `GSE125105` | 699 | 2796 |
| `GSE130051` | 1501 | 1501 |
| `GSE132203` | 795 | 1590 |
| `GSE134379` | 808 | 2424 |
| `GSE140686` | 1504 | 1504 |
| `GSE141065` | 557 | 1114 |
| `GSE145361` | 1889 | 4708 |
| `GSE147221` | 720 | 2139 |
| `GSE147740` | 1128 | 3384 |
| `GSE153712` | 726 | 1923 |
| `GSE157131` | 1218 | 2436 |
| `GSE179325` | 574 | 1148 |
| `GSE183647` | 565 | 1130 |
| `GSE183920` | 603 | 2412 |
| `GSE185920` | 1471 | 2942 |
| `GSE197678` | 2922 | 5844 |
| `GSE208713` | 568 | 1664 |
| `GSE210255` | 1394 | 1394 |
| `GSE215240` | 937 | 0 |
| `GSE219037` | 553 | 1755 |
| `GSE223817` | 984 | 2307 |
| `GSE224124` | 1107 | 1107 |
| `GSE224365` | 712 | 2130 |
| `GSE225845` | 595 | 1050 |
| `GSE226569` | 694 | 1388 |
| `GSE270375` | 994 | 994 |
| `GSE280206` | 648 | 1296 |
| `GSE40279` | 138 | 276 |
| `GSE42861` | 241 | 810 |
| `GSE43414` | 50 | 61 |
| `GSE51032` | 297 | 594 |
| `GSE55763` | 1841 | 5523 |
| `GSE56046` | 290 | 290 |
| `GSE56105` | 168 | 504 |
| `GSE59685` | 11 | 22 |
| `GSE65362` | 562 | 562 |
| `GSE68379` | 1028 | 0 |
| `GSE72680` | 12 | 36 |
| `GSE73801` | 585 | 1093 |
| `GSE74193` | 55 | 30 |
| `GSE77716` | 573 | 572 |
| `GSE78874` | 16 | 32 |
| `GSE80417` | 675 | 1986 |

## Census delta

- Unique GSM before (disk census): **157941**
- Unique GSM after: **157993**
- GEO phenotype rows before merge (in-memory): **0**
- GEO observed rows in previous disk census: **33961**
- GEO observed phenotype rows after: **74971**
- Authoritative phenotype rows added this merge: **74971**

> **Dirty disk baseline:** previous `census.json` already had GEO rows. That does **not** mean this merge added zero — use `n_phenotype_rows_added` / in-memory before_merge (0 on a full rebuild).

GEO backfill does not add `sample` rows (EWAS_db scan does). Unique-GSM movement is EWAS_db mirror growth, not GEO. Phenotype-row movement is the GEO delta.

## Operator notes

- GEO rows are omitted for Hub GSM (Hub wins).
- Atlas enrichment stays on study.metadata_json only.
- Disease/cancer rows need explicit case/control tokens.
- Diagnosis-only text stays in metadata_json.
- Training heads read Hub pack Parquet, not geo_metadata_backfill.
- Authoritative GEO Δ is merge_stats.n_phenotype_rows_added (in-memory phenotypes before merge are usually 0 on a full refresh).
- For a clean incremental test: MBS_SKIP_GEO_BACKFILL=1 refresh → assert zero geo_metadata_backfill rows → fetch/merge → compare exact Δ.
