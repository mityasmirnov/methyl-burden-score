# Raw data snapshot (2026-08-05)

Sanitized inspection of remote EWAS DataHub hosts and example files under
`$MBS_DATA_ROOT/raw`. Machine-readable details:
[`summary.json`](summary.json).

## Remote hosts

| Role | NGDC-advertised FTP | Working HTTP mirror (this server) |
|------|---------------------|-----------------------------------|
| **All Data** | `ftp://download.big.ac.cn/ewas/datahub/EWAS_db/` | `https://download.cncb.ac.cn/ewas/datahub/EWAS_db/` |
| **Baseline Data** | `ftp://download.big.ac.cn/ewas/datahub/download/` | `https://download.cncb.ac.cn/ewas/datahub/download/` |

Notes:

- Direct FTP from `power-horse` stalls on connect; prefer the HTTP mirrors or FileZilla.
- **All Data** index lists **~1989** study directories (GSE / E-MTAB / ENCODE / CPTAC / …), dated 13-Apr-2026.
- **Baseline Data** HTTP packs match the NGDC table (`*_v1.zip` names), plus `GMQN.zip`.

### All Data layout (`EWAS_db/{STUDY}/`)

Per-study directories contain one `GSM*.txt` per sample. Example
(`GSE104451/GSM2800663.txt`):

```text
cg00050873	0.852
cg00212031	0.046
...
```

Two-column TSV: Illumina probe ID, beta value. No header.

Sample counts (HTTP listing):

| Study | `GSM*.txt` files |
|-------|------------------|
| GSE104451 | 61 |
| GSE42861 | 689 |
| GSE87571 | 732 |
| GSE35069 | 60 |

### Stage 0 labeling GSEs in All Data

**All 18** requested Stage 0 GSEs are present under `EWAS_db/` (including those
absent from public CpGCorpus S3: e.g. GSE42861, GSE87571, GSE104451, GSE108423).

### Baseline Data packs (HTTP sizes sampled)

| File | Content-Length |
|------|----------------|
| `tissue_methylation_v1.zip` | 8,266,723,841 (~7.7 GiB) |
| `sample_tissue_methylation_v1.zip` | 195,615 |
| `brain_methylation_v1.zip` | 2,975,575,853 |
| `blood_methylation_v1.zip` | 5,220,215,339 |

Full pack list: [`docs/EWAS_DATA.md`](../../EWAS_DATA.md).

## Local `$MBS_DATA_ROOT/raw` (at inspection time)

Approximate total **~7.9 GiB** on disk.

| Path | Size | Status |
|------|------|--------|
| `cpgcorpus/` | ~6.7 GiB | Stage 0 GSEs complete; plus leftovers from aborted full sync (~46 other GSE dirs) |
| `ewas_atlas/` | ~269 MiB | **Complete** (5 batch files) |
| `ewas_datahub/` | ~482 MiB | `tissue_methylation_v1.zip` partial (~6% of 7.7 GiB); `EWAS_db/` empty (FTP stall) |
| `manifests/epicv2/` | ~446 MiB | **Complete** EPICv2 reannotated v3.0 csv.gz |

### EWAS Atlas examples

| File | Bytes | Schema (head) |
|------|-------|---------------|
| `EWAS_Atlas_associations.tsv` | 106,146,695 | `Association_ID`, `probe_ID`, `trait`, `case_beta`, `control_beta`, `correlation`, `p_value`, `study_ID`, `PMID` (13 cols) |
| `EWAS_Atlas_studies.tsv` | ~151 KiB | `study_ID`, `trait`, `case_description`, `control_description`, `PMID` |
| `EWAS_Atlas_cohorts.tsv` | ~389 KiB | `cohort_ID`, `study_ID`, `platform`, `sample_size`, ages, `tissue`, `ancestry`, … |
| `EWAS_Atlas_probe_annotations.tsv` | ~166 MiB | `Probe.id`, `Chr`, `Pos`, `Related.transcripts`, `Pos.to.TSS`, `CpG.islands`, … |
| `EWAS_trait_trait_logP.txt` | ~964 KiB | trait×trait matrix (first column = trait name) |

### CpGCorpus Stage 0 examples (Arrow)

Orientation: **samples as rows**, probes as columns, plus `GSM_ID`.

| File | Rows | Cols | Notes |
|------|------|------|-------|
| `GSE35069/.../gse_betas.arrow` | 60 | 485,578 | GEO betas (450K-scale) |
| `GSE35069/.../metadata.arrow` | 60 | 34 | GEO sample metadata strings |
| `GSE116992/.../QCDPB.arrow` | 22 | 865,920 | sesame QCDPB betas (EPIC) |
| `GSE116992/.../metadata.arrow` | 22 | 38 | |
| `GSE125367/.../QCDPB.arrow` | 44 | 865,920 | |
| `GSE125367/.../metadata.arrow` | 44 | 36 | |

### EPICv2 manifest

`EPICv2_reannotated_manifest_v3.0.csv.gz` (~446 MiB): Illumina-style columns
starting `IlmnID`, `Name`, probe sequences, design fields, plus GENCODEv49 /
GeneHancer / clock annotations (Zenodo v3.0).

### DataHub baseline zip (local)

`tissue_methylation_v1.zip` was still downloading at inspection; partial archive
is not yet zip-readable.

## Implications

1. Prefer **HTTP** mirrors for DataHub All Data / Baseline on this host.
2. Stage 0 GSEs missing from CpGCorpus S3 can be pulled from **EWAS_db** as
   per-sample beta text files (different format than CpGCorpus Arrow).
3. Keep Atlas (associations) and DataHub (matrices) as separate raw trees.

## Download method update (2026-08-05)

HTTP root https://download.cncb.ac.cn/ewas/datahub/ is mirrored with an HTML-index parser (`scripts/download_ewas_datahub.sh`) because recursive wget drops `GSM*.txt` under JS-enhanced listings. Trees: `EWAS_db/`, `add_ewas_db/`, `download/`.
