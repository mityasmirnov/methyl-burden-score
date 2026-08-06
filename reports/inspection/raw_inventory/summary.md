# Raw data inventory (for systematic organization)

Inspected: **2026-08-06**. Machine-readable details:
[`summary.json`](summary.json).

Purpose: make it obvious **what landed where**, **what format each tree uses**,
and **how to reorganize** without mixing sources.

## Totals under `data/raw/`

| Tree | Bytes (approx.) | Role |
|------|-----------------|------|
| `ewas_datahub/` | **~74 GiB** | DataHub All Data + baseline packs (in progress) |
| `cpgcorpus/` | **~6.7 GiB** | CpGCorpus Arrow (Stage 0 + sync leftovers) |
| `manifests/` | **~446 MiB** | EPICv2 reannotated manifest |
| `ewas_atlas/` | **~269 MiB** | Atlas association knowledge (complete) |
| **all raw** | **~81 GiB** | |

## Recommended canonical layout

```text
data/raw/
  cpgcorpus/{GSE}/{GPL}/{betas,metadata}/*.arrow     # training matrices (Arrow)
  ewas_atlas/*.tsv                                   # association knowledge
  ewas_datahub/
    EWAS_db/{STUDY}/*.txt                            # per-sample probe β text
    download/*_methylation_v1.zip                    # baseline trait bundles
    download/sample_*.zip                            # sample metadata (often R)
    download/GMQN.zip
    add_ewas_db/                                     # supplemental (currently empty)
  manifests/epicv2/*.csv.gz                          # array annotation only
```

Do **not** mix Atlas TSVs, DataHub zips, and CpGCorpus Arrow into one folder.
Tag provenance in the catalog (`source_system`: `cpgcorpus` | `ewas_atlas` |
`ewas_datahub_db` | `ewas_datahub_baseline` | `epicv2_manifest`).

---

## 1. EWAS Atlas — `data/raw/ewas_atlas/` (complete)

| File | Bytes | Data rows | Cols | Schema (head) |
|------|------:|----------:|-----:|---------------|
| `EWAS_Atlas_associations.tsv` | 106,146,695 | **804,919** | 13 | `Association_ID`, `probe_ID`, `trait`, `case_description`, `case_beta`, `control_description`, `control_beta`, `correlation`, `p_value`, `rank_in_study`, `effect_size`, `study_ID`, `PMID` |
| `EWAS_Atlas_studies.tsv` | 153,681 | **1,902** | 5 | `study_ID`, `trait`, `case_description`, `control_description`, `PMID` |
| `EWAS_Atlas_cohorts.tsv` | 397,573 | **3,983** | 15 | `cohort_ID`, `study_ID`, `stage`, `platform`, `sample_size`, ages, `tissue`, `ancestry`, … |
| `EWAS_Atlas_probe_annotations.tsv` | 174,062,150 | **900,413** | 8 | `Probe.id`, `Chr`, `Pos`, `Related.transcripts`, `Pos.to.TSS`, `Related.GWAS.sites`, `Pos.to.GWAS.sites`, `CpG.islands` |
| `EWAS_trait_trait_logP.txt` | 986,305 | **371** (+header) | ~373 | trait×trait matrix; first column = trait name |

**Top association row (truncated):**

```text
EA00332384  cg21197425  fetal intolerance of labor  …  pos  …
```

**Top study row:**

```text
ES00033  body mass index (BMI)  …  PMID 24630777
```

**Top cohort row:**

```text
52  ES00033  Initial  450K  sample_size=239  male_percentage=0.85  …
```

**Top probe row:**

```text
cg00000029  chr16  53468112  RBL2(...)  TSS offsets  …
```

---

## 2. EPICv2 manifest — `data/raw/manifests/epicv2/` (complete)

| File | Bytes | Data rows | Cols |
|------|------:|----------:|-----:|
| `EPICv2_reannotated_manifest_v3.0.csv.gz` | 466,768,830 | **937,691** | **62** |

Columns start Illumina-style (`IlmnID`, `Name`, probe sequences, design fields)
then GENCODEv49 / GeneHancer / clock annotations (Zenodo v3.0).

**Top row (truncated):**

```text
cg25324105_BC11  cg25324105  1754126  ATTTATAAACTAATAACCCAAAATACATTTCCCAAAAACCTTCACAACCA  …
```

---

## 3. CpGCorpus — `data/raw/cpgcorpus/`

### Stage 0 selected (keep)

Orientation: **samples as rows**, probes as columns, key `GSM_ID`.

| GSE | Bytes | Platform | Betas file | Samples×probes | Metadata |
|-----|------:|----------|------------|----------------:|----------|
| GSE116992 | 422 MiB | GPL13534 | `QCDPB.arrow` | 7 × 485,579 | 7 × 38 |
| GSE116992 | (same) | GPL21145 | `QCDPB.arrow` | 22 × 865,920 | 22 × 38 |
| GSE125367 | 415 MiB | GPL21145 | `QCDPB.arrow` | 44 × 865,920 | 44 × 36 |
| GSE35069 | 313 MiB | GPL13534 | `gse_betas.arrow` | 60 × 485,578 | 60 × 34 |

Metadata columns include GEO fields: `GSM_ID`, `title`, `geo_accession`,
`source_name_ch1`, `characteristics_ch1`, …

### Leftovers (reorganize)

**46 extra GSE directories** remain from the aborted full-corpus sync
(e.g. GSE100184…GSE110184). They are **not** the Stage 0 labeling set.

Suggested cleanup:

```text
data/raw/cpgcorpus/          # keep only Stage 0 GSEs
data/raw/cpgcorpus/_partial_fullsync/   # move the 46 leftovers here
```

---

## 4. EWAS DataHub — `data/raw/ewas_datahub/`

HTTP root: https://download.cncb.ac.cn/ewas/datahub/

### 4a. Baseline packs — `download/` (partial)

| File | Bytes on disk | Zip readable? | Members | Notes |
|------|--------------:|:-------------:|--------:|-------|
| `GMQN.zip` | 42,008,375 | yes | 16 | complete |
| `age_methylation_v1.zip` | 12,598,215,716 | yes | 2 | present |
| `ancestry_category_methylation_v1.zip` | 2,109,764,152 | yes | 2 | present |
| `blood_methylation_v1.zip` | 5,220,215,339 | yes | 2 | complete vs HTTP length |
| `bmi_methylation_v1.zip` | 3,284,635,193 | yes | 2 | present |
| `brain_methylation_v1.zip` | 2,975,575,853 | yes | 2 | complete |
| `cancer_methylation_v1.zip` | ~0.5–GB+ | **no** | — | still downloading |

**Organized 2026-08-06:** unique flat packs moved into `download/`;
identical flats (blood/brain) deleted; incomplete flat age deleted
(download/ copy larger). No `ewas_datahub/*.zip` remain.

### 4b. Flat duplicates under `ewas_datahub/*.zip` — done

Previously flats sat beside `download/`. Consolidated:

| Flat file | Action |
|-----------|--------|
| `tissue_methylation_v1.zip`, `sex_methylation_v1.zip`, `sample_*.zip` | moved → `download/` |
| `blood_methylation_v1.zip`, `brain_methylation_v1.zip` | deleted (byte-identical to `download/`) |
| `age_methylation_v1.zip` | deleted incomplete flat; kept larger `download/` copy |

### 4c. Sample metadata (unzips)

Remote `sample_*_methylation_v1.zip` members include **`RDX3`** RData plus a
space-quoted `.txt` table. Prefer `.txt`.

**Current local layout (2026-08-06):** all nine families unpacked under
`reports/inspection/ewas_datahub_samples/` (zips deleted after extract).
Ancestry members are named `sample_race.*`. Structure profiles:
`reports/inspection/ewas_metadata_structure/` and [`docs/EWAS_METADATA.md`](../../docs/EWAS_METADATA.md).

### 4d. All Data — `EWAS_db/` (in progress)

| Metric | Value |
|--------|------:|
| Study dirs on disk | 359 (many **empty** leftovers from failed recursive wget) |
| Non-empty studies | **2** so far (`CGCI-HTMCP-CC`, `CPTAC-3`) |
| Sample files | **1,665** |
| Bytes | **~22 GiB** |

Remote has ~1989 studies; download job is early (parser walking study-by-study).

**File format (every sample):** no header, two columns `probe_id<TAB>beta`:

```text
cg18478105	0.04
cg09835024	0.094
cg14361672	0.659
…
```

**Naming is heterogeneous:**

- GEO-like: `GSM….txt`
- CPTAC: UUID-like or `C3L-….txt` / `C3N-….txt` (sometimes `_935k` twin)
- CGCI: `HTMCP-….txt`

Junk artifact `(.+?)` under `CGCI-HTMCP-CC/` **deleted** (2026-08-06).
Empty leftover study dirs from failed recursive wget **removed** (357);
non-empty remain: `CGCI-HTMCP-CC`, `CPTAC-3` (All Data pull still running).

### 4e. `add_ewas_db/`

Directories `add_txt_450/`, `add_txt_850/`, `add_txt_935/` exist but were
**empty** on the HTTP server at pull time (+ `move.sh` only).

### 4f. CpGCorpus quarantine

Stage 0 kept under `raw/cpgcorpus/`: `GSE116992`, `GSE125367`, `GSE35069`.
**46** leftover GSE dirs from the aborted full sync moved to
`raw/cpgcorpus/_partial_fullsync/` (see README there). Not for default ingest.

---

## Organization checklist

1. ~~**Consolidate DataHub baseline** into `ewas_datahub/download/` only~~ **done**
2. ~~**Quarantine CpGCorpus leftovers** (`_partial_fullsync/`) vs Stage 0 GSEs~~ **done**
3. ~~**Delete empty `EWAS_db` dirs** and the `(.+?)` junk file~~ **done**; resume All Data pull until Stage 0 GSEs are present as text betas (or keep using CpGCorpus Arrow where available).
4. ~~**Catalog provenance lanes**~~ **done** — `provenance_lane` + `source_release.source_system` (`cpgcorpus` | `ewas_atlas` | `ewas_datahub_db` | `ewas_datahub_baseline` | `epicv2_manifest`). Assay-level `file_format` / `sha256` remain on `assay_file` as ingest fills them.
5. ~~**Keep Atlas vs DataHub vs CpGCorpus separate**~~ **done** (lanes + raw roots).

## Related docs

- [`docs/EWAS_DATA.md`](../../docs/EWAS_DATA.md) — download policy & hosts  
- [`docs/CPGCORPUS_STAGE0.md`](../../docs/CPGCORPUS_STAGE0.md) — Stage 0 GSE list  
- [`reports/inspection/raw_data_snapshot/summary.md`](../raw_data_snapshot/summary.md) — earlier host/layout notes  
