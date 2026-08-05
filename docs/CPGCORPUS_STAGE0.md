# CpGCorpus Stage 0 GSE selection

Selective downloads go to `$MBS_DATA_ROOT/raw/cpgcorpus/{GSE}/…` via:

```bash
source scripts/activate_data_environment.sh
export AWS_DEFAULT_REGION=us-east-1
bash scripts/download_cpgcorpus_gse.sh
# or: bash scripts/download_cpgcorpus_gse.sh GSE116992 GSE35069
# or: make download-cpgcorpus-gse
```

Full-corpus sync remains `scripts/download_cpgcorpus.sh` (large; use only when needed).

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

Availability was checked against `s3://cpgpt-lucascamillo-public/data/cpgcorpus/raw/` (requester-pays, `us-east-1`). Studies marked **no** are not present under that prefix and must be obtained from GEO/other archives if needed for Stage 0.

## Layout (per GSE)

```text
data/raw/cpgcorpus/{GSE}/{GPL}/
  betas/QCDPB.arrow          # or gse_betas.arrow
  metadata/metadata.arrow
```

## Provenance

- Bucket: `s3://cpgpt-lucascamillo-public/data/cpgcorpus/raw`
- GSE list: `configs/data/stage0_cpgcorpus_gse_list.txt`
- Downloader: `scripts/download_cpgcorpus_gse.sh`
