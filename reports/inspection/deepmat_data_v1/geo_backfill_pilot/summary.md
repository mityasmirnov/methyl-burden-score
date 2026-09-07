# GEO backfill pilot summary

- Generated: `2026-09-07T15:04:43Z`
- GEO parquet GSM in: **56201**
- Catalog samples touched (EWAS_db-only): **48394**
- Hub-skipped GSM (already complete): **6985**
- Hub GSM with GEO fill-missing: **92**
- Species-quarantined GSM (non-human/unknown): **0**
- GEO GSM not in catalog: **822**
- Phenotype rows added: **98636**
- Samples with ≥1 observed GEO phenotype: **43853**
- Samples with `metadata_json.geo`: **48394**
- Multi-study GSM (membership): **0**

## Invariants

- GEO phenotype rows on Hub pack-membership GSM: **0** (must be 0)
- Atlas blobs on `sample.metadata_json`: **0** (must be 0)
- `sample_source_membership` includes `geo_metadata_backfill`: **False** (must be false)

## Tissue ontology

- mapped / unmapped / ambiguous / empty: **43263** / **12116** / **0** / **0**
- Unmapped examples: `LM`, `PR`, `sarcoma`, `colon cancer cell line`, `ovarian cancer cell line`, `embryonic stem cell, undifferentiated`, `embryonic stem cell, differentiated;spontaneous.day3`, `embryonic stem cell, differentiated;NeuralProgenitors.day22`, `embryonic stem cell, differentiated;OligodendrocyteProgenitors.day42`, `parthenogentic embryonic stem cell, undifferentiated`, `induced pluripotent stem cell, undifferentiated`, `induced pluripotent stem cell, differentiated;NeuralProgenitors.day22`, `induced pluripotent stem cell, differentiated;NeuralProgenitors.day23`, `induced pluripotent stem cell, differentiated;OligodendrocyteProgenitors.day42`, `Mole`, `Somatic.Primary`, `Somatic.Tissue`, `brain_whole`, `brain_neg`, `brain_pos`

## Phenotypes by id

| phenotype_id | rows | observed | unique GSM |
| --- | ---: | ---: | ---: |
| `age` | 21940 | 21940 | 21940 |
| `cancer` | 1687 | 1687 | 1687 |
| `disease` | 6065 | 6065 | 6065 |
| `sex` | 32083 | 32083 | 32083 |
| `tissue` | 36861 | 36861 | 36861 |

## Label status

| phenotype_id | label_status | n |
| --- | --- | ---: |
| `age` | `observed` | 21940 |
| `cancer` | `case` | 1422 |
| `cancer` | `control` | 265 |
| `disease` | `case` | 1937 |
| `disease` | `control` | 4128 |
| `sex` | `observed` | 32083 |
| `tissue` | `observed` | 36861 |

## Eligibility (`source_family=geo_metadata_backfill`)

| phenotype_id | n | cases | controls | core | aux |
| --- | ---: | ---: | ---: | --- | --- |
| `age` | 21940 | None | None | True | True |
| `cancer` | 1687 | 1422 | 265 | True | True |
| `disease` | 6065 | 1937 | 4128 | True | True |
| `sex` | 32083 | None | None | False | True |
  - `sex` not core: sex is auxiliary biological / QC, not a core burden target
| `tissue` | 36861 | None | None | True | True |

## Per study (merge)

| study_id | samples touched | phenotype rows |
| --- | ---: | ---: |
| `GSE103286` | 112 | 392 |
| `GSE103332` | 116 | 396 |
| `GSE104210` | 675 | 1920 |
| `GSE105018` | 1230 | 1230 |
| `GSE107351` | 116 | 116 |
| `GSE107352` | 51 | 51 |
| `GSE108187` | 102 | 0 |
| `GSE108576` | 89 | 141 |
| `GSE109379` | 1104 | 1104 |
| `GSE110184` | 70 | 0 |
| `GSE110544` | 58 | 58 |
| `GSE110697` | 78 | 78 |
| `GSE111165` | 147 | 289 |
| `GSE111629` | 316 | 1060 |
| `GSE114598` | 61 | 0 |
| `GSE115278` | 474 | 1422 |
| `GSE116339` | 679 | 2037 |
| `GSE117859` | 608 | 1824 |
| `GSE119144` | 60 | 60 |
| `GSE122038` | 217 | 203 |
| `GSE122126` | 101 | 165 |
| `GSE122920` | 104 | 104 |
| `GSE122994` | 380 | 760 |
| `GSE123601` | 142 | 111 |
| `GSE124617` | 94 | 5 |
| `GSE125105` | 699 | 2796 |
| `GSE125895` | 269 | 793 |
| `GSE130051` | 1501 | 1501 |
| `GSE131013` | 240 | 720 |
| `GSE131482` | 153 | 219 |
| `GSE132203` | 795 | 2385 |
| `GSE132650` | 78 | 78 |
| `GSE132804` | 334 | 668 |
| `GSE134379` | 808 | 2424 |
| `GSE135205` | 51 | 0 |
| `GSE136790` | 108 | 108 |
| `GSE137845` | 54 | 0 |
| `GSE140686` | 1504 | 0 |
| `GSE141039` | 153 | 282 |
| `GSE141065` | 557 | 1114 |
| `GSE141254` | 84 | 202 |
| `GSE142512` | 395 | 1185 |
| `GSE143708` | 75 | 150 |
| `GSE144977` | 89 | 89 |
| `GSE145361` | 1889 | 4708 |
| `GSE147221` | 720 | 2852 |
| `GSE147430` | 132 | 0 |
| `GSE147740` | 1128 | 4512 |
| `GSE153712` | 726 | 1923 |
| `GSE157131` | 1218 | 2436 |
| `GSE157397` | 51 | 49 |
| `GSE179325` | 574 | 1148 |
| `GSE183647` | 565 | 1130 |
| `GSE183920` | 603 | 2412 |
| `GSE185920` | 1471 | 2942 |
| `GSE197678` | 2922 | 8766 |
| `GSE208713` | 568 | 1096 |
| `GSE210255` | 1394 | 2788 |
| `GSE215240` | 937 | 937 |
| `GSE219037` | 553 | 1755 |
| `GSE223817` | 984 | 2307 |
| `GSE224124` | 1107 | 1107 |
| `GSE224218` | 170 | 340 |
| `GSE224365` | 712 | 2130 |
| `GSE225845` | 595 | 1645 |
| `GSE226569` | 694 | 1388 |
| `GSE237561` | 126 | 378 |
| `GSE240482` | 58 | 0 |
| `GSE267015` | 68 | 0 |
| `GSE270375` | 994 | 1988 |
| `GSE280206` | 648 | 1296 |
| `GSE31848` | 153 | 153 |
| `GSE34387` | 76 | 0 |
| `GSE40279` | 138 | 276 |
| `GSE40699` | 62 | 0 |
| `GSE41336` | 90 | 0 |
| `GSE42861` | 241 | 810 |
| `GSE43414` | 142 | 61 |
| `GSE49177` | 61 | 51 |
| `GSE51032` | 297 | 891 |
| `GSE52025` | 62 | 62 |
| `GSE55763` | 1841 | 5523 |
| `GSE56046` | 290 | 580 |
| `GSE56105` | 168 | 504 |
| `GSE57204` | 54 | 0 |
| `GSE57342` | 338 | 210 |
| `GSE57853` | 92 | 125 |
| `GSE59685` | 11 | 22 |
| `GSE65183` | 144 | 0 |
| `GSE65362` | 562 | 562 |
| `GSE66872` | 60 | 0 |
| `GSE68060` | 120 | 0 |
| `GSE68379` | 1028 | 0 |
| `GSE68825` | 144 | 533 |
| `GSE70783` | 74 | 149 |
| `GSE72680` | 12 | 36 |
| `GSE73801` | 585 | 1416 |
| `GSE74193` | 55 | 85 |
| `GSE74738` | 79 | 63 |
| `GSE77348` | 85 | 0 |
| `GSE77716` | 573 | 1145 |
| `GSE78874` | 16 | 48 |
| `GSE79185` | 61 | 12 |
| `GSE80417` | 675 | 1986 |
| `GSE81472` | 51 | 0 |
| `GSE86078` | 149 | 149 |
| `GSE86355` | 75 | 0 |
| `GSE87177` | 119 | 86 |
| `GSE90496` | 2801 | 2801 |
| `GSE94943` | 72 | 24 |

## Census delta

- Unique GSM before (disk census): **173076**
- Unique GSM after: **173076**
- GEO phenotype rows before merge (in-memory): **0**
- GEO observed rows in previous disk census: **110077**
- GEO observed phenotype rows after: **98636**
- Authoritative phenotype rows added this merge: **98636**

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
