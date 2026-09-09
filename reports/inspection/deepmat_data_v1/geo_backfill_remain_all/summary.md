# Next GEO backfill GSE list

- Generated: `2026-09-08T15:02:40Z`
- Candidates (≥1 GSM, not fetched): **1345**
- Selected: **1345** → `configs/data/geo_backfill_remain_all_gse.txt`

| study_id | n_assay | soft | age | tissue_mapped | score |
| --- | ---: | --- | ---: | ---: | ---: |
| `GSE51057` | 329 | True | 329 | 0 | 987 |
| `GSE137898` | 265 | True | 265 | 0 | 795 |
| `GSE137903` | 254 | True | 254 | 0 | 762 |
| `GSE137841` | 250 | True | 135 | 0 | 520 |
| `GSE137682` | 152 | False | None | None | 152 |
| `GSE113018` | 138 | False | None | None | 138 |
| `GSE156968` | 138 | False | None | None | 138 |
| `GSE79144` | 138 | False | None | None | 138 |
| `GSE79257` | 137 | False | None | None | 137 |
| `GSE89251` | 136 | False | None | None | 136 |
| `GSE71955` | 135 | False | None | None | 135 |
| `GSE113967` | 134 | False | None | None | 134 |
| `GSE157252` | 134 | False | None | None | 134 |
| `GSE39672` | 133 | False | None | None | 133 |
| `GSE104293` | 132 | False | None | None | 132 |
| `GSE54503` | 132 | False | None | None | 132 |
| `GSE88929` | 132 | False | None | None | 132 |
| `GSE159526` | 131 | False | None | None | 131 |
| `GSE245924` | 131 | False | None | None | 131 |
| `GSE168779` | 130 | False | None | None | 130 |
| `GSE226257` | 130 | False | None | None | 130 |
| `GSE235414` | 130 | False | None | None | 130 |
| `GSE59457` | 130 | False | None | None | 130 |
| `GSE155333` | 129 | False | None | None | 129 |
| `GSE63704` | 129 | False | None | None | 129 |
| `GSE128601` | 128 | False | None | None | 128 |
| `GSE240184` | 128 | False | None | None | 128 |
| `GSE107459` | 127 | False | None | None | 127 |
| `GSE226823` | 127 | False | None | None | 127 |
| `GSE89181` | 127 | False | None | None | 127 |
| `GSE52731` | 126 | False | None | None | 126 |
| `GSE93933` | 126 | False | None | None | 126 |
| `GSE79009` | 125 | False | None | None | 125 |
| `GSE228100` | 124 | False | None | None | 124 |
| `GSE189778` | 123 | False | None | None | 123 |
| `GSE149747` | 122 | False | None | None | 122 |
| `GSE218549` | 122 | False | None | None | 122 |
| `GSE74432` | 122 | False | None | None | 122 |
| `GSE87095` | 122 | False | None | None | 122 |
| `GSE101961` | 121 | False | None | None | 121 |

## Notes

- Species gate: series with any non-human GSM (when SOFT cached) are omitted.
- Uncached SOFTs are ranked by assay N only; species checked at fetch time.
- Age/tissue richness is unknown until SOFT fetch for uncached studies.
- Larger crawl remains gated for training; this only picks the next fetch list.
