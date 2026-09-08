# GEO backfill pilot summary

- Generated: `2026-09-07T18:10:23Z`
- GEO parquet GSM in: **92601**
- Catalog samples touched (EWAS_db-only): **80035**
- Hub-skipped GSM (already complete): **10361**
- Hub GSM with GEO fill-missing: **2159**
- Species-quarantined GSM (non-human/unknown): **357**
- GEO GSM not in catalog: **1848**
- Phenotype rows added: **149653**
- Samples with ≥1 observed GEO phenotype: **70756**
- Samples with `metadata_json.geo`: **80035**
- Multi-study GSM (membership): **329**

## Invariants

- GEO phenotype rows on Hub pack-membership GSM: **208** (must be 0)
- Atlas blobs on `sample.metadata_json`: **0** (must be 0)
- `sample_source_membership` includes `geo_metadata_backfill`: **False** (must be false)

## Tissue ontology

- mapped / unmapped / ambiguous / empty: **59877** / **30519** / **0** / **0**
- Unmapped examples: `2621`, `2622`, `2691`, `2692`, `3391`, `3392`, `2791`, `2792`, `3161`, `3162`, `3211`, `3212`, `3431`, `3432`, `3441`, `3442`, `3311`, `3312`, `3341`, `3342`

## Phenotypes by id

| phenotype_id | rows | observed | unique GSM |
| --- | ---: | ---: | ---: |
| `age` | 32969 | 32969 | 32969 |
| `cancer` | 3438 | 3438 | 3438 |
| `disease` | 10516 | 10516 | 10516 |
| `sex` | 52413 | 52413 | 52413 |
| `tissue` | 50317 | 50317 | 50317 |

## Label status

| phenotype_id | label_status | n |
| --- | --- | ---: |
| `age` | `observed` | 32969 |
| `cancer` | `case` | 3015 |
| `cancer` | `control` | 423 |
| `disease` | `case` | 4315 |
| `disease` | `control` | 6201 |
| `sex` | `observed` | 52413 |
| `tissue` | `observed` | 50317 |

## Eligibility (`source_family=geo_metadata_backfill`)

| phenotype_id | n | cases | controls | core | aux |
| --- | ---: | ---: | ---: | --- | --- |
| `age` | 32969 | None | None | True | True |
| `cancer` | 3438 | 3015 | 423 | True | True |
| `disease` | 10516 | 4315 | 6201 | True | True |
| `sex` | 52413 | None | None | False | True |
  - `sex` not core: sex is auxiliary biological / QC, not a core burden target
| `tissue` | 50317 | None | None | True | True |

## Per study (merge)

| study_id | samples touched | phenotype rows |
| --- | ---: | ---: |
| `GSE100227` | 479 | 880 |
| `GSE100264` | 386 | 1158 |
| `GSE101764` | 112 | 336 |
| `GSE103253` | 245 | 52 |
| `GSE104210` | 675 | 1920 |
| `GSE105018` | 1230 | 1230 |
| `GSE105109` | 384 | 880 |
| `GSE106648` | 89 | 356 |
| `GSE107080` | 405 | 1215 |
| `GSE107298` | 394 | 0 |
| `GSE107351` | 116 | 116 |
| `GSE107352` | 51 | 51 |
| `GSE108187` | 102 | 0 |
| `GSE108576` | 89 | 141 |
| `GSE109379` | 1104 | 1104 |
| `GSE110128` | 233 | 92 |
| `GSE110184` | 70 | 0 |
| `GSE110544` | 58 | 58 |
| `GSE110697` | 78 | 78 |
| `GSE111165` | 147 | 289 |
| `GSE111629` | 316 | 1060 |
| `GSE112611` | 402 | 1280 |
| `GSE112893` | 374 | 748 |
| `GSE113904` | 232 | 464 |
| `GSE114598` | 61 | 0 |
| `GSE115278` | 474 | 1422 |
| `GSE116339` | 679 | 2037 |
| `GSE117859` | 608 | 1824 |
| `GSE117860` | 529 | 1058 |
| `GSE119144` | 60 | 60 |
| `GSE121633` | 480 | 1440 |
| `GSE122038` | 217 | 203 |
| `GSE122126` | 101 | 165 |
| `GSE122920` | 104 | 104 |
| `GSE122994` | 380 | 760 |
| `GSE123601` | 142 | 111 |
| `GSE124413` | 500 | 959 |
| `GSE124617` | 94 | 5 |
| `GSE125105` | 699 | 2796 |
| `GSE125895` | 269 | 793 |
| `GSE128821` | 536 | 536 |
| `GSE130051` | 1501 | 1501 |
| `GSE131013` | 240 | 480 |
| `GSE131433` | 440 | 440 |
| `GSE131482` | 153 | 219 |
| `GSE131989` | 371 | 371 |
| `GSE132181` | 391 | 782 |
| `GSE132203` | 795 | 2385 |
| `GSE132650` | 78 | 78 |
| `GSE132804` | 334 | 668 |
| `GSE134379` | 808 | 2424 |
| `GSE135205` | 51 | 0 |
| `GSE136790` | 108 | 108 |
| `GSE137841` | 250 | 385 |
| `GSE137845` | 54 | 0 |
| `GSE137898` | 265 | 530 |
| `GSE137903` | 254 | 508 |
| `GSE138279` | 65 | 195 |
| `GSE140686` | 1504 | 0 |
| `GSE141039` | 153 | 282 |
| `GSE141065` | 557 | 1114 |
| `GSE141254` | 84 | 168 |
| `GSE142512` | 395 | 1185 |
| `GSE144129` | 454 | 454 |
| `GSE144858` | 300 | 996 |
| `GSE144977` | 89 | 89 |
| `GSE145361` | 1889 | 4708 |
| `GSE146376` | 280 | 560 |
| `GSE147221` | 720 | 2852 |
| `GSE147430` | 132 | 0 |
| `GSE147740` | 1128 | 4512 |
| `GSE149412` | 292 | 584 |
| `GSE150643` | 240 | 480 |
| `GSE151042` | 492 | 484 |
| `GSE151732` | 256 | 512 |
| `GSE153712` | 726 | 1923 |
| `GSE154566` | 233 | 466 |
| `GSE156984` | 244 | 244 |
| `GSE157131` | 1218 | 2436 |
| `GSE157341` | 274 | 548 |
| `GSE157397` | 51 | 49 |
| `GSE163970` | 492 | 984 |
| `GSE174422` | 256 | 768 |
| `GSE175458` | 547 | 1641 |
| `GSE179325` | 574 | 1148 |
| `GSE183647` | 565 | 1130 |
| `GSE183920` | 603 | 2412 |
| `GSE185008` | 230 | 230 |
| `GSE185061` | 369 | 738 |
| `GSE185920` | 1471 | 2942 |
| `GSE197512` | 372 | 0 |
| `GSE197678` | 2922 | 8766 |
| `GSE197723` | 317 | 647 |
| `GSE201322` | 232 | 696 |
| `GSE207460` | 290 | 870 |
| `GSE208529` | 391 | 391 |
| `GSE208713` | 568 | 1096 |
| `GSE210254` | 418 | 836 |
| `GSE210255` | 1394 | 2788 |
| `GSE212682` | 320 | 640 |
| `GSE215240` | 937 | 937 |
| `GSE219037` | 553 | 1755 |
| `GSE220622` | 391 | 782 |
| `GSE221029` | 230 | 0 |
| `GSE221219` | 320 | 0 |
| `GSE221256` | 305 | 609 |
| `GSE223467` | 321 | 321 |
| `GSE223817` | 984 | 2307 |
| `GSE224124` | 1107 | 1107 |
| `GSE224218` | 170 | 340 |
| `GSE224339` | 239 | 80 |
| `GSE224365` | 712 | 2130 |
| `GSE225810` | 241 | 74 |
| `GSE225845` | 595 | 1645 |
| `GSE226569` | 694 | 1388 |
| `GSE232332` | 275 | 693 |
| `GSE234461` | 245 | 735 |
| `GSE237103` | 380 | 380 |
| `GSE237561` | 126 | 378 |
| `GSE240482` | 58 | 0 |
| `GSE240704` | 363 | 363 |
| `GSE243075` | 363 | 253 |
| `GSE243529` | 524 | 524 |
| `GSE246337` | 500 | 1500 |
| `GSE249113` | 398 | 1194 |
| `GSE250513` | 376 | 736 |
| `GSE252169` | 359 | 718 |
| `GSE253101` | 411 | 573 |
| `GSE263434` | 285 | 0 |
| `GSE263787` | 298 | 298 |
| `GSE267015` | 68 | 0 |
| `GSE269244` | 238 | 238 |
| `GSE269416` | 290 | 289 |
| `GSE270223` | 392 | 392 |
| `GSE270375` | 994 | 1988 |
| `GSE280206` | 648 | 1296 |
| `GSE31848` | 153 | 153 |
| `GSE34387` | 76 | 0 |
| `GSE40279` | 138 | 276 |
| `GSE40699` | 62 | 0 |
| `GSE41336` | 90 | 0 |
| `GSE42861` | 241 | 810 |
| `GSE43414` | 142 | 61 |
| `GSE49177` | 9 | 9 |
| `GSE51032` | 297 | 891 |
| `GSE52025` | 62 | 62 |
| `GSE52270` | 269 | 321 |
| `GSE52401` | 244 | 112 |
| `GSE53740` | 187 | 530 |
| `GSE55763` | 1841 | 5523 |
| `GSE56046` | 290 | 580 |
| `GSE56105` | 168 | 504 |
| `GSE56588` | 243 | 382 |
| `GSE57204` | 54 | 0 |
| `GSE57342` | 338 | 210 |
| `GSE57853` | 92 | 125 |
| `GSE59250` | 434 | 85 |
| `GSE59685` | 11 | 22 |
| `GSE60185` | 285 | 730 |
| `GSE64509` | 4 | 8 |
| `GSE65183` | 144 | 0 |
| `GSE65362` | 562 | 562 |
| `GSE66872` | 60 | 0 |
| `GSE68060` | 120 | 0 |
| `GSE68379` | 1028 | 0 |
| `GSE68825` | 144 | 533 |
| `GSE69550` | 268 | 0 |
| `GSE70783` | 74 | 149 |
| `GSE71678` | 343 | 0 |
| `GSE72680` | 12 | 36 |
| `GSE72773` | 148 | 444 |
| `GSE72774` | 317 | 1062 |
| `GSE72775` | 171 | 513 |
| `GSE73103` | 203 | 609 |
| `GSE73801` | 585 | 1416 |
| `GSE74104` | 267 | 180 |
| `GSE74193` | 55 | 85 |
| `GSE74738` | 79 | 63 |
| `GSE75248` | 335 | 1 |
| `GSE77348` | 85 | 0 |
| `GSE77696` | 378 | 378 |
| `GSE77716` | 573 | 1145 |
| `GSE78874` | 16 | 48 |
| `GSE79185` | 61 | 12 |
| `GSE80417` | 675 | 1986 |
| `GSE80970` | 8 | 16 |
| `GSE81472` | 51 | 0 |
| `GSE82273` | 730 | 1460 |
| `GSE84207` | 330 | 330 |
| `GSE84493` | 310 | 0 |
| `GSE84727` | 609 | 974 |
| `GSE85210` | 253 | 235 |
| `GSE85212` | 763 | 763 |
| `GSE86078` | 149 | 149 |
| `GSE86355` | 75 | 0 |
| `GSE87177` | 119 | 86 |
| `GSE87571` | 433 | 1727 |
| `GSE87640` | 17 | 17 |
| `GSE87648` | 69 | 207 |
| `GSE89353` | 620 | 1032 |
| `GSE90496` | 2801 | 2801 |
| `GSE93646` | 428 | 428 |
| `GSE94943` | 72 | 24 |
| `GSE97362` | 91 | 352 |
| `GSE97628` | 380 | 189 |

## Census delta

- Unique GSM before (disk census): **173076**
- Unique GSM after: **173076**
- GEO phenotype rows before merge (in-memory): **0**
- GEO observed rows in previous disk census: **97382**
- GEO observed phenotype rows after: **149653**
- Authoritative phenotype rows added this merge: **149653**

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
