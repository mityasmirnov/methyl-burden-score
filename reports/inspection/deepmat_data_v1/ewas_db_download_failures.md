# EWAS_db download failures

- Generated: `2026-09-02T10:41:52Z`
- Log: `/data/projects/methyl-burden-score/artifacts/logs/downloads/ewas_datahub_EWAS_db.log`
- EWAS_db root: `/data/projects/methyl-burden-score/data/raw/ewas_datahub/EWAS_db`
- Studies with ≥1 logged failure: **1475**
- Total `WARN: failed` lines: **5136**
- Still missing or empty on disk: **4224**
- HTML-parse artifact filenames (`(.+?)`): **1474**
- Last study progress in log: **1582** / **1989**

## Retry

Manifest: `artifacts/logs/downloads/ewas_db_retry_manifest.tsv`

```bash
bash scripts/retry_ewas_db_download_failures.sh
```

## Top studies by failure count

| study_id | failures | still_missing |
| --- | ---: | ---: |
| `GSE197678` | 254 | 200 |
| `GSE112611` | 95 | 62 |
| `GSE56046` | 81 | 61 |
| `GSE140686` | 79 | 63 |
| `GSE59685` | 59 | 52 |
| `GSE141441` | 57 | 55 |
| `GSE210255` | 57 | 39 |
| `CPTAC-3` | 52 | 39 |
| `GSE163970` | 51 | 42 |
| `GSE174422` | 50 | 38 |
| `GSE157131` | 48 | 32 |
| `GSE55763` | 48 | 30 |
| `GSE68379` | 43 | 38 |
| `GSE168779` | 42 | 31 |
| `GSE185920` | 37 | 32 |
| `GSE112893` | 36 | 24 |
| `GSE59250` | 32 | 26 |
| `GSE183647` | 31 | 29 |
| `GSE51032` | 31 | 26 |
| `GSE130051` | 30 | 23 |
| `GSE169646` | 30 | 23 |
| `GSE210254` | 30 | 24 |
| `GSE141065` | 28 | 18 |
| `GSE153712` | 27 | 25 |
| `GSE161476` | 24 | 20 |

## Notes

- Re-run `bash scripts/download_ewas_datahub.sh EWAS_db` to resume; successful files are skipped via `wget -c`.
- Post-download hook runs `mbs catalog refresh-release` automatically (disable with `EWAS_DATAHUB_SKIP_POST_HOOK=1`).
