# CpGCorpus Stage 0 GSE selection

Nothing in this repository auto-downloads CpGCorpus. Downloads run only when a
script below is invoked explicitly. All files land under
`$MBS_DATA_ROOT/raw/cpgcorpus`; logs under `$MBS_ARTIFACT_ROOT/logs/downloads`.

## Foreground downloads

Selective Stage 0 GSE sync:

```bash
source scripts/activate_data_environment.sh
export AWS_DEFAULT_REGION=us-east-1
bash scripts/download_cpgcorpus_gse.sh
# or: bash scripts/download_cpgcorpus_gse.sh GSE116992 GSE35069
# or: make download-cpgcorpus-gse
```

Full-corpus sync remains `scripts/download_cpgcorpus.sh` (large; use only when needed).

## Background downloads (nohup)

```bash
source scripts/activate_data_environment.sh
bash scripts/download_cpgcorpus_background.sh gse
# or: bash scripts/download_cpgcorpus_background.sh gse GSE116992
# or: bash scripts/download_cpgcorpus_background.sh full
```

The wrapper prints a PID and log path. Monitor with `tail -f` on the nohup log
under `$MBS_ARTIFACT_ROOT/logs/downloads/`.

## Requested labeling cohort

| GSE | Labeling expectation | In public CpGCorpus S3? |
|-----|----------------------|-------------------------|
| GSE116992 | BAFopathy / CSS + NCBRS | yes (`GPL13534`, `GPL21145`) |
| GSE66552 | 7q11.23 del/dup; Williams + Dup7 | **no** |
| GSE74432 | Sotos / NSD1 | **no** |
| GSE97362 | CHARGE + Kabuki | **no** |
| GSE116300 | Kabuki phenotype | **no** |
| GSE95040 | ICF subtypes + controls | **no** |
| GSE104451 | Silver-Russell | **no** |
| GSE125367 | NCBRS / SMARCA2 | yes (`GPL21145`) |
| GSE55491 | Silver-Russell | **no** |
| GSE108423 | KDM5C LOF / ID (CJS in AJHG table) | **no** |
| GSE89353 | Unsolved congenital disorders / NDD | **no** |
| GSE52588 | Down syndrome blood; DS vs mother/sibling | **no** |
| GSE42861 | Rheumatoid arthritis cases/controls | **no** |
| GSE85210 | Healthy smokers vs never smokers | **no** |
| GSE87571 | Aging cohort, whole blood | **no** |
| GSE87648 | IBD cases/controls | **no** |
| GSE99863 | ENID trial Gambian children | **no** |
| GSE35069 | Purified blood cell types (FlowSorted-style) | yes (`GPL13534`) |

Availability was checked against `s3://cpgpt-lucascamillo-public/data/cpgcorpus/raw/` (requester-pays, `us-east-1`). Studies marked **no** are not present under that prefix.

**Alternate source:** all 18 GSEs above **are** present in EWAS DataHub All Data
(`https://download.cncb.ac.cn/ewas/datahub/EWAS_db/{GSE}/` as per-sample
`GSM*.txt` beta files). See [`EWAS_DATA.md`](EWAS_DATA.md) and
[`reports/inspection/raw_data_snapshot/summary.md`](../reports/inspection/raw_data_snapshot/summary.md).

## Layout (per GSE)

```text
data/raw/cpgcorpus/{GSE}/{GPL}/
  betas/QCDPB.arrow          # or gse_betas.arrow
  metadata/metadata.arrow
```

## Provenance

- Bucket: `s3://cpgpt-lucascamillo-public/data/cpgcorpus/raw`
- GSE list: `configs/data/stage0_cpgcorpus_gse_list.txt`
- Downloaders: `scripts/download_cpgcorpus_gse.sh`, `scripts/download_cpgcorpus.sh`
- Background wrapper: `scripts/download_cpgcorpus_background.sh`
