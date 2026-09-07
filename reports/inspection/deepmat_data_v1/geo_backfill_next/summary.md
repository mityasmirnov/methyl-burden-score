# Next GEO backfill GSE list

- Generated: `2026-09-07T14:51:24Z`
- Candidates (≥50 GSM, not fetched): **539**
- Selected: **100** → `configs/data/geo_backfill_next_gse.txt`

| study_id | n_assay | soft | age | tissue_mapped | score |
| --- | ---: | --- | ---: | ---: | ---: |
| `GSE82273` | 889 | False | None | None | 889 |
| `GSE84727` | 847 | False | None | None | 847 |
| `GSE85212` | 763 | False | None | None | 763 |
| `GSE87571` | 732 | False | None | None | 732 |
| `GSE89353` | 620 | False | None | None | 620 |
| `GSE175458` | 547 | False | None | None | 547 |
| `GSE128821` | 536 | False | None | None | 536 |
| `GSE117860` | 529 | False | None | None | 529 |
| `GSE243529` | 524 | False | None | None | 524 |
| `GSE72774` | 508 | False | None | None | 508 |
| `GSE124413` | 500 | False | None | None | 500 |
| `GSE246337` | 500 | False | None | None | 500 |
| `GSE151042` | 492 | False | None | None | 492 |
| `GSE163970` | 492 | False | None | None | 492 |
| `GSE121633` | 480 | False | None | None | 480 |
| `GSE100227` | 479 | False | None | None | 479 |
| `GSE144129` | 454 | False | None | None | 454 |
| `GSE131433` | 440 | False | None | None | 440 |
| `GSE138279` | 439 | False | None | None | 439 |
| `GSE59250` | 434 | False | None | None | 434 |
| `GSE93646` | 428 | False | None | None | 428 |
| `GSE210254` | 418 | False | None | None | 418 |
| `GSE253101` | 411 | False | None | None | 411 |
| `GSE107080` | 405 | False | None | None | 405 |
| `GSE112611` | 402 | False | None | None | 402 |
| `GSE249113` | 398 | False | None | None | 398 |
| `GSE107298` | 394 | False | None | None | 394 |
| `GSE270223` | 392 | False | None | None | 392 |
| `GSE132181` | 391 | False | None | None | 391 |
| `GSE208529` | 391 | False | None | None | 391 |
| `GSE220622` | 391 | False | None | None | 391 |
| `GSE100264` | 386 | False | None | None | 386 |
| `GSE105109` | 384 | False | None | None | 384 |
| `GSE53740` | 384 | False | None | None | 384 |
| `GSE87648` | 384 | False | None | None | 384 |
| `GSE237103` | 380 | False | None | None | 380 |
| `GSE97628` | 380 | False | None | None | 380 |
| `GSE77696` | 378 | False | None | None | 378 |
| `GSE250513` | 376 | False | None | None | 376 |
| `GSE112893` | 374 | False | None | None | 374 |

## Notes

- Species gate: series with any non-human GSM (when SOFT cached) are omitted.
- Uncached SOFTs are ranked by assay N only; species checked at fetch time.
- Age/tissue richness is unknown until SOFT fetch for uncached studies.
- Larger crawl remains gated for training; this only picks the next fetch list.
