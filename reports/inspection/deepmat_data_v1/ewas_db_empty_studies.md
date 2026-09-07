# EWAS_db dirs with no GSM*.txt

- Generated: `2026-09-07T08:58:06Z`
- Count: **339** / 1989 study dirs
- Truly empty or artifact-only: **294** (artifact-only: 123)
- Dirs that already have non-GSM sample `.txt`: **45**

## By accession family

- `ArrayExpress`: 6
- `GSE`: 292
- `TCGA`: 33
- `consortium_other`: 7
- `other`: 1

## Cause

mirror_ewas_db only wget files matching ^GSM[0-9]+\.txt$. TCGA/ArrayExpress/ENCODE/CPTAC use other sample filenames; those dirs that already have files came from an earlier broader mirror. Empty GSE dirs are failed/missing GSM downloads or HTML-index parse failures (artifact (.+?)).

## Remote probes (sample)

| study_id | remote GSM | remote other txt | sample |
| --- | ---: | ---: | --- |
| `CGCI-HTMCP-CC` | 0 | 120 | HTMCP-03-06-02001-01A.txt, HTMCP-03-06-02002-01A.txt, HTMCP-03-06-02003-01A.txt |
| `CPTAC-3` | 0 | 1998 | 2f2e5477-42a4-4906-a943-bf7f80.txt, 5a84eae1-197e-4463-ad65-59becc.txt, 93e30fd5-e57e-4503-a175-863c7d.txt |
| `E-MTAB-1274` | 0 | 20 | CB11130.txt, CB11171.txt, CB152.txt |
| `E-MTAB-1866` | 0 | 648 | TWPID10016.txt, TWPID10017.txt, TWPID10024.txt |
| `GSE102119` | 146 | 0 | GSM2724241.txt, GSM2724242.txt, GSM2724243.txt |
| `GSE107143` | 16 | 0 | GSM2861662.txt, GSM2861663.txt, GSM2861664.txt |
| `GSE110724` | 21 | 0 | GSM3015190.txt, GSM3015191.txt, GSM3015192.txt |
| `GSE84399` | 0 | 1 | GSE84399.txt |
| `TARGET-ALL-P3` | 0 | 112 | TARGET-10-PANJRD-09A.txt, TARGET-10-PANKMB-09A.txt, TARGET-10-PANRIM-03A.txt |
| `TARGET-AML` | 0 | 197 | TARGET-20-PABLDZ-14A-01D.txt, TARGET-20-PADYIR-04A-01D.txt, TARGET-20-PADZCG-04A-02D.txt |
| `TCGA-ACC` | 0 | 80 | TCGA-OR-A5J1-01A-11D-A29J-05.txt, TCGA-OR-A5J2-01A-11D-A29J-05.txt, TCGA-OR-A5J3-01A-11D-A29J-05.txt |
| `TCGA-BLCA` | 0 | 440 | TCGA-2F-A9KO-01A-11D-A38H-05.txt, TCGA-2F-A9KP-01A-11D-A38H-05.txt, TCGA-2F-A9KQ-01A-11D-A38H-05.txt |

## Follow-ups

1. Keep GSM-preferring filter for GSE studies (matches Hub GSM IDs).
2. `mirror_ewas_db` now falls back to non-GSM `*.txt` when a study has no GSM
   (TCGA / ArrayExpress / ENCODE / CPTAC). Re-run mirror or a targeted pass for
   those dirs.
3. Re-run `make retry-ewas-db-failures` for empty GSE dirs that still have remote GSM
   (e.g. GSE102119 has 146 remote GSM but local dir empty).
4. `scan_ewas_db_tree` / convert skip parse artifacts (`(.+?)`); non-GSM samples on
   disk are inventoried. Delete leftover `(.+?)` files opportunistically.
