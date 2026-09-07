# EWAS_db download failures

- Generated: `2026-09-07T14:23:14Z`
- Log: `/data/projects/methyl-burden-score/artifacts/logs/downloads/ewas_datahub_EWAS_db.log`
- EWAS_db root: `/data/projects/methyl-burden-score/data/raw/ewas_datahub/EWAS_db`
- Studies with ≥1 logged failure: **1818**
- Total `WARN: failed` lines: **8142**
- Still missing or empty on disk: **0**
- HTML-parse artifact filenames (`(.+?)`): **1818**
- Last study progress in log: **1989** / **1989**

## Retry

Manifest: `artifacts/logs/downloads/ewas_db_retry_manifest.tsv`

```bash
bash scripts/retry_ewas_db_download_failures.sh
```

## Top studies by failure count

| study_id | failures | still_missing |
| --- | ---: | ---: |
| `GSE90496` | 274 | 0 |
| `GSE197678` | 254 | 0 |
| `GSE89353` | 110 | 0 |
| `GSE112611` | 95 | 0 |
| `GSE84727` | 85 | 0 |
| `GSE82273` | 82 | 0 |
| `GSE56046` | 81 | 0 |
| `GSE87571` | 81 | 0 |
| `GSE140686` | 79 | 0 |
| `GSE87648` | 61 | 0 |
| `GSE59685` | 59 | 0 |
| `GSE73801` | 58 | 0 |
| `GSE141441` | 57 | 0 |
| `GSE210255` | 57 | 0 |
| `GSE74193` | 56 | 0 |
| `CPTAC-3` | 52 | 0 |
| `GSE163970` | 51 | 0 |
| `GSE174422` | 50 | 0 |
| `GSE77716` | 49 | 0 |
| `TCGA-PRAD` | 49 | 0 |
| `GSE157131` | 48 | 0 |
| `GSE55763` | 48 | 0 |
| `GSE72774` | 45 | 0 |
| `TCGA-BLCA` | 45 | 0 |
| `GSE85212` | 44 | 0 |

## Notes

- Re-run `bash scripts/download_ewas_datahub.sh EWAS_db` to resume; successful files are skipped via `wget -c`.
- Post-download hook runs `make catalog-refresh-release` automatically (Atlas GSE map seed + census; disable with `EWAS_DATAHUB_SKIP_POST_HOOK=1`).
