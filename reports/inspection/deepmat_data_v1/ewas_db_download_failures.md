# EWAS_db download failures

- Generated: `2026-09-02T11:38:44Z`
- Log: `/data/projects/methyl-burden-score/artifacts/logs/downloads/ewas_datahub_EWAS_db.log`
- EWAS_db root: `/data/projects/methyl-burden-score/data/raw/ewas_datahub/EWAS_db`
- Studies with ≥1 logged failure: **1475**
- Total `WARN: failed` lines: **5138**
- Still missing or empty on disk: **2710**
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
| `GSE197678` | 254 | 199 |
| `GSE112611` | 95 | 61 |
| `GSE56046` | 81 | 60 |
| `GSE140686` | 79 | 62 |
| `GSE59685` | 59 | 51 |
| `GSE141441` | 57 | 54 |
| `GSE210255` | 57 | 38 |
| `CPTAC-3` | 52 | 0 |
| `GSE163970` | 51 | 41 |
| `GSE174422` | 50 | 37 |
| `GSE157131` | 48 | 31 |
| `GSE55763` | 48 | 29 |
| `GSE68379` | 43 | 37 |
| `GSE168779` | 42 | 30 |
| `GSE185920` | 37 | 31 |
| `GSE112893` | 36 | 23 |
| `GSE59250` | 32 | 25 |
| `GSE183647` | 31 | 28 |
| `GSE51032` | 31 | 25 |
| `GSE130051` | 30 | 22 |
| `GSE169646` | 30 | 22 |
| `GSE210254` | 30 | 23 |
| `GSE141065` | 28 | 17 |
| `GSE153712` | 27 | 24 |
| `GSE161476` | 24 | 19 |

## Notes

- Re-run `bash scripts/download_ewas_datahub.sh EWAS_db` to resume; successful files are skipped via `wget -c`.
- Post-download hook runs `mbs catalog refresh-release` automatically (disable with `EWAS_DATAHUB_SKIP_POST_HOOK=1`).
