# Milestone 7B convert progress

Updated: `2026-08-24T16:09:28Z`

**7B packs:** 5 / 6 done

| Family | Status | Phase | Samples | Phenotype rows | Loci | betas.zarr |
|--------|--------|-------|--------:|---------------:|-----:|-----------:|
| `ancestry` | `done` | `done` | 1380 | 1380 | 482387 | 1.2 GiB |
| `bmi` | `done` | `done` | 2070 | 2070 | 482387 | 1.8 GiB |
| `brain` | `done` | `done` | 1997 | 1997 | 482387 | 1.7 GiB |
| `blood` | `done` | `done` | 3402 | 3402 | 482387 | 2.9 GiB |
| `cancer` | `done` | `done` | 10101 | 10841 | 482387 | 9.1 GiB |
| `disease` | `pending` | `pending` | — | — | — | — |

- Running convert-pack: `disease`
- Virtual index present: `False`
- Inspection summary present: `False`

Phases: `pending` → `stream_scratch` → `write_zarr` → `qc_or_checksum` → `done`.

## Remaining

- Finish convert for: disease (skip-if-exists resume)
- mbs matrix index-hub-packs --check-overlap
- scripts/write_stage0_7b_report.py
- mbs catalog refresh-release (7A pointers only)
- Required checks; mark TODO_PIPELINE 7B done with evidence
