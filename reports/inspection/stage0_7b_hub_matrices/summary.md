# Milestone 7B: complete Hub pack matrices

Per-pack Zarr stores (no dense nine-pack union). Age/tissue/sex full
matrices were not reconverted.

## 7B matrices

| Family | Matrix ID | Samples | Loci | Phenotype rows | Platform | Content sha256 |
|--------|-----------|--------:|-----:|---------------:|----------|----------------|
| `disease` | `matrix-hub-disease-full-v1` | 12218 | 482387 | 14501 | `450K` | `57a8ea341c08…` |
| `cancer` | `matrix-hub-cancer-full-v1` | 10101 | 482387 | 10841 | `450K` | `3427719c63fa…` |
| `blood` | `matrix-hub-blood-full-v1` | 3402 | 482387 | 3402 | `450K` | `630523780b2d…` |
| `brain` | `matrix-hub-brain-full-v1` | 1997 | 482387 | 1997 | `450K` | `dd7ef36c8e15…` |
| `bmi` | `matrix-hub-bmi-full-v1` | 2070 | 482387 | 2070 | `450K` | `7fc617293599…` |
| `ancestry` | `matrix-hub-ancestry-full-v1` | 1380 | 482387 | 1380 | `450K` | `e7c78230840b…` |

Disease/cancer: `n_samples` is unique GSM (matrix rows);
`n_phenotype_rows` is long-form sample-info (may exceed unique GSM).

## Virtual index

- Path: `canonical/matrices/hub_pack_matrix_index.parquet`
- Present: True
- Rows / unique GSM / families: 47843 / 34234 / 9

## Overlap concordance

- Status: `concordant`
- Shared GSM: 9409
- Pairs checked / discordant: 18346 / 0
- Max abs diff: 0.0
- Merge allowed: True

## Frozen 5d matrices (not overwritten)

- `matrix-hub-age-full-v1` shape `[8374, 482379]` platform `HM450`
- `matrix-hub-tissue-full-v1` shape `[5323, 482379]` platform `HM450`
- `matrix-hub-sex-full-v1` shape `[2978, 482379]` platform `HM450`
