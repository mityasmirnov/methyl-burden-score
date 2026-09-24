# Orphan RBS qualification census

CPU census only — does **not** change assignment. Product Stage B / OOF should keep one column per **qualified multi-CpG** orphan `region_id`; singletons belong on the direct path.

Graph: `graph-grch38-gencode38-cgi-tile-v2`
max_loci: `65536`
gene_allocation: `explicit_only`

| Matrix | n_orphan | singleton (1 CpG) | multi (≥2 CpG) | n_direct | n_genes |
|--------|---------:|------------------:|---------------:|---------:|--------:|
| `matrix-hub-age-tissue-sex-full-v1` | 1953 | 861 | 1092 | 9322 | 2646 |
| `matrix-hub-nine-pack-virtual-v1` | 1953 | 861 | 1092 | 9322 | 2646 |

## Cardinality detail

### `matrix-hub-age-tissue-sex-full-v1`

| CpGs/region | n_orphan_regions |
|------------:|-----------------:|
| 1 | 861 |
| 2 | 378 |
| 3 | 287 |
| 4 | 188 |
| 5 | 110 |
| 6 | 53 |
| 7 | 32 |
| 8 | 17 |
| 9 | 10 |
| 10 | 2 |
| 11 | 1 |
| 12 | 3 |
| 13 | 1 |
| 14 | 2 |
| 15 | 1 |
| 16 | 1 |
| 19 | 1 |
| 20 | 1 |
| 22 | 1 |
| 23 | 1 |
| 27 | 1 |
| 41 | 1 |

### `matrix-hub-nine-pack-virtual-v1`

| CpGs/region | n_orphan_regions |
|------------:|-----------------:|
| 1 | 861 |
| 2 | 378 |
| 3 | 287 |
| 4 | 188 |
| 5 | 110 |
| 6 | 53 |
| 7 | 32 |
| 8 | 17 |
| 9 | 10 |
| 10 | 2 |
| 11 | 1 |
| 12 | 3 |
| 13 | 1 |
| 14 | 2 |
| 15 | 1 |
| 16 | 1 |
| 19 | 1 |
| 20 | 1 |
| 22 | 1 |
| 23 | 1 |
| 27 | 1 |
| 41 | 1 |

