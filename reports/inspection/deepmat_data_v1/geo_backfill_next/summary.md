# Next GEO backfill GSE list

- Generated: `2026-09-08T09:39:58Z`
- Candidates (≥50 GSM, not fetched): **440**
- Selected: **100** → `configs/data/geo_backfill_next_gse.txt`

| study_id | n_assay | soft | age | tissue_mapped | score |
| --- | ---: | --- | ---: | ---: | ---: |
| `GSE51057` | 329 | True | 329 | 0 | 987 |
| `GSE161476` | 229 | False | None | None | 229 |
| `GSE58218` | 228 | False | None | None | 228 |
| `GSE56600` | 227 | False | None | None | 227 |
| `GSE279073` | 226 | False | None | None | 226 |
| `GSE147040` | 221 | False | None | None | 221 |
| `GSE72021` | 221 | False | None | None | 221 |
| `GSE214568` | 220 | False | None | None | 220 |
| `GSE53051` | 220 | False | None | None | 220 |
| `GSE156994` | 219 | False | None | None | 219 |
| `GSE124366` | 215 | False | None | None | 215 |
| `GSE185090` | 215 | False | None | None | 215 |
| `GSE56581` | 214 | False | None | None | 214 |
| `GSE104942` | 210 | False | None | None | 210 |
| `GSE114134` | 205 | False | None | None | 205 |
| `GSE229203` | 200 | False | None | None | 200 |
| `GSE60753` | 199 | False | None | None | 199 |
| `GSE220838` | 196 | False | None | None | 196 |
| `GSE281307` | 196 | False | None | None | 196 |
| `GSE49149` | 196 | False | None | None | 196 |
| `GSE113725` | 194 | False | None | None | 194 |
| `GSE40576` | 194 | False | None | None | 194 |
| `GSE269983` | 192 | False | None | None | 192 |
| `GSE36054` | 192 | False | None | None | 192 |
| `GSE60132` | 192 | False | None | None | 192 |
| `GSE77718` | 192 | False | None | None | 192 |
| `GSE103186` | 191 | False | None | None | 191 |
| `GSE66351` | 190 | False | None | None | 190 |
| `GSE230850` | 189 | False | None | None | 189 |
| `GSE75067` | 188 | False | None | None | 188 |
| `GSE199591` | 187 | False | None | None | 187 |
| `GSE143752` | 186 | False | None | None | 186 |
| `GSE69138` | 185 | False | None | None | 185 |
| `GSE69270` | 184 | False | None | None | 184 |
| `GSE66836` | 183 | False | None | None | 183 |
| `GSE227809` | 182 | False | None | None | 182 |
| `GSE256394` | 182 | False | None | None | 182 |
| `GSE103659` | 181 | False | None | None | 181 |
| `GSE108462` | 181 | False | None | None | 181 |
| `GSE143755` | 181 | False | None | None | 181 |

## Notes

- Species gate: series with any non-human GSM (when SOFT cached) are omitted.
- Uncached SOFTs are ranked by assay N only; species checked at fetch time.
- Age/tissue richness is unknown until SOFT fetch for uncached studies.
- Larger crawl remains gated for training; this only picks the next fetch list.
