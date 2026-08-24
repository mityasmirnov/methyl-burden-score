# Data-source inspection guide

## Principle

No source enters a canonical release because its filename or documentation looks plausible. Every GSE/GPL, DataHub bundle, annotation file, vocabulary, and model feature artifact is inspected, checksummed, and summarized before conversion.

Coding agents should inspect sanitized reports under `reports/inspection/`, not raw sample-level data.

## Standard source report

Each source produces:

```text
reports/inspection/{source_id}/
├── summary.md
├── summary.json
├── schema.json
├── file_inventory.parquet
├── sample_alignment.parquet
├── value_qc.parquet
├── metadata_counts.parquet
├── mapping_qc.parquet
├── graph_coverage.parquet
└── warnings.json
```

These files contain schemas and aggregate statistics. Do not write identifiable or sensitive sample fields into the report.

## File inventory

For every source file, record:

```text
source accession
absolute path under /data
format
byte size
SHA-256
retrieval timestamp
source URL or command
processing level
genome build
expected matrix role
```

Example shell inventory:

```bash
find "$MBS_DATA_ROOT/raw" -type f -printf '%p\t%s\n' | sort
sha256sum /data/path/to/source-file
```

## CpGCorpus and Arrow files

CpGCorpus is typically organized by GSE and GPL with beta and metadata Arrow files. Inspect every combination independently.

For the Stage 0 labeling cohort and which GSEs exist on the public requester-pays bucket, see [`CPGCORPUS_STAGE0.md`](CPGCORPUS_STAGE0.md). Prefer `scripts/download_cpgcorpus_gse.sh` over a full-corpus sync unless the entire corpus is required.

Minimal Python inspection:

```python
from pathlib import Path

import pyarrow.feather as feather

path = Path("/data/path/to/QCDPB.arrow")
table = feather.read_table(path, memory_map=True)
print(table.schema)
print(table.num_rows, table.num_columns)
print(table.column_names[:20])
```

Determine rather than assume:

- whether samples are rows or columns;
- how sample IDs are encoded;
- whether feature names are Illumina probes or genomic coordinates;
- whether metadata are embedded in the beta table;
- whether duplicate field names exist;
- whether Arrow list columns or dense numeric columns are used.

## Sample and metadata alignment

For each GSE/GPL report:

```text
number of beta samples
number of metadata samples
intersection size
beta-only IDs
metadata-only IDs
duplicate sample IDs
duplicate donor IDs
technical replicate groups
longitudinal groups
```

Stop conversion when unexplained rows remain unmatched.

## Beta-value checks

Report globally and per sample:

- minimum, maximum, mean, standard deviation, and quantiles;
- fraction outside `[0,1]`;
- missing-value fraction;
- zero-variance loci;
- duplicate loci and duplicate probe IDs;
- identical or near-identical sample vectors;
- extreme sample means;
- expected versus observed platform probe count.

Do not silently clip out-of-range values or replace missing values during inspection.

## Age inspection

Record:

```text
numeric conversion rate
missingness
minimum and maximum
median and quantiles
likely units
prenatal or neonatal encodings
implausible values
age distribution by study, tissue, sex, platform, and processing level
```

Age transformations are fitted only inside training folds.

## Tissue inspection

Record:

```text
raw label frequencies
normalized ontology mapping
ambiguous mappings
classes represented by one study only
classes represented by one platform only
classes with fewer than the configured minimum samples
```

A class confined to one study cannot support a claim of study-independent tissue prediction.

## Confounding tables

Always produce aggregate tables for:

```text
trait × study
trait × platform
trait × processing level
trait × sex
trait × age band
study × platform
study × tissue
```

Also train or evaluate coverage-only and platform-only controls before interpreting model performance.

## Probe-to-locus mapping

Report:

- probes mapped to one canonical GRCh38 locus;
- unmapped probes;
- ambiguous probe mappings;
- probes with genome-build disagreement;
- multiple probes collapsed to one locus;
- loci represented on 450K, EPIC, and EPICv2;
- loci absent from each platform.

Coordinate convention and cytosine position must be explicit.

## Annotation graph inspection

For each graph release report:

```text
loci, regions, genes, and edges
regions by biological role
loci per region and per gene
multi-gene locus edges
unassigned/intergenic loci
coverage by platform and study
coverage by region role
MethylGPT vocabulary overlap
coordinate or contig mismatches
```

Manually inspect representative genes with simple and complex transcript structures before freezing the graph.

## CpGPT static-feature inspection

Verify:

```text
source repository and commit
model and checkpoint hash
DNA model and context length
raw and adapted dimensions
coordinate-to-row mapping
NaN and infinity count
vector norm distribution
per-dimension variance
sampled equality with direct encode_sequence calls
```

The exporter must not invoke sample encoding or dynamic transformer states.

## MethylGPT vocabulary and token-prior inspection

Verify:

```text
vocabulary checksum
special-token order
number of probe tokens
duplicate probe IDs
probe-to-GRCh38 mapping rate
ambiguous mappings
coverage by platform, gene, and region
embedding-table shape
checkpoint compatibility
```

A missing or incompatible checkpoint key is an error; partial loading is not acceptable for canonical export.

## EWAS DataHub

Keep DataHub external validation / Hub-pack training releases separate from
CpGCorpus until study overlap and normalization provenance are resolved.

**Policy: download all public DataHub data** (All Data `EWAS_db` plus Baseline
packs). See [`EWAS_DATA.md`](EWAS_DATA.md).

Inspect:

- study and GSE identifiers (`project_id` in sample-info; study dirs under `EWAS_db`);
- trait bundle definitions (baseline `*_methylation_v1.zip` packs);
- **sample information** (see below);
- matrix processing and batch-correction description (GMQN);
- platform and genome build;
- overlap with CpGCorpus studies;
- phenotype coding and primary columns per family.

### Sample-info packs (required before Hub phenotype joins / Milestone 5c)

Remote archives are `sample_*_methylation_v1.zip` under
`raw/ewas_datahub/download/`. After download they are **unzipped once** to
Cursor-visible extracts:

```text
reports/inspection/ewas_datahub_samples/sample_{family}_methylation_v1/
  sample_{family}.txt     # R write.table (prefer this)
  sample_{family}.RData   # binary sidecar (do not index)
```

**Exception:** ancestry zip is `sample_ancestry_category_methylation_v1.zip` but
members are `sample_race.txt` / `sample_race.RData`.

Nine families are profiled: age, ancestry, blood, bmi, brain, cancer, disease,
sex, tissue. Contracts (parse recipe, join keys, family→column map, overlaps):

- [`EWAS_METADATA.md`](EWAS_METADATA.md)
- `reports/inspection/ewas_metadata_structure/` (`mbs inspect ewas-metadata`)

```bash
uv run mbs inspect ewas-metadata
make export-ewas-sample-info FAMILY=tissue
# → $MBS_DATA_ROOT/canonical/phenotypes/{family}_sample_info.parquet
```

Do not invent Hub column names; use `FAMILY_VALUE_COLUMN` in
`mbs.registry.sample_info`. Each zip also contains RData (`RDX3`); prefer the
`.txt` member.

Shallow download layout / sizes / sample Ns:
[`reports/inspection/raw_inventory/`](../reports/inspection/raw_inventory/) and
[`DATA_CATALOG.md`](DATA_CATALOG.md) (refresh with
`uv run python scripts/write_raw_inventory_refresh.py`).
## Model performance reports (Stage 0)

Sanitized train/eval summaries (no sample-level betas):

| Report | Milestone |
|--------|-----------|
| `reports/inspection/stage0_5b_benchmark/` | 5b fixture holdouts |
| `reports/inspection/stage0_hub_real_benchmark/` | 5b″ Hub study-holdout packs |
| `reports/inspection/stage0_5d_max_n/` | 5d max-N age/tissue/sex DeepRVAT flat |

Regenerate the 5d report after a finished run:

```bash
uv run python scripts/write_stage0_5d_report.py
```

Raw `metrics.jsonl` / checkpoints stay under `$MBS_ARTIFACT_ROOT` and are not
committed.

## EWAS Atlas

Atlas is an association evidence source, not a sample matrix. Freeze and
checksum association, study, cohort, and annotation exports. Map CpGs through
the project graph rather than relying only on supplied nearest-gene fields.
Full batch exports: [`EWAS_DATA.md`](EWAS_DATA.md).

**Small tables (Cursor-indexed; inspect before 5c):**

| File | Role |
|------|------|
| `EWAS_Atlas_studies.tsv` | study_ID (ES…), trait, PMID |
| `EWAS_Atlas_cohorts.tsv` | cohort_ID → study_ID, platform, ages, tissue, ancestry |
| `EWAS_trait_trait_logP.txt` | trait×trait logP matrix |

Encoding is latin-1; rare malformed cohort rows with extra tabs are skipped by
`mbs inspect ewas-metadata`. **Do not** join Atlas `study_ID` to Hub
`project_id` by string equality (ES* vs GSE*).

Large tables (associations, probe annotations) stay out of the Cursor index;
profile them only via sanitized reports when needed. Structure contracts:
[`EWAS_METADATA.md`](EWAS_METADATA.md).

## Long-read and methylartist inputs

For BAM or methylartist SQLite sources verify:

```text
BAM sorted and indexed
reference contigs agree with GRCh38
MM/ML tags and modification codes
primary and supplementary alignment fractions
mapping-quality distribution
duplicate calls at one read/locus
methylated, unmethylated, and no-call counts
coverage and strand balance
SQLite table schema and genomic indexes
```

Long-read calls are reduced to the same canonical locus IDs, with coverage and count fields retained.

## Acceptance gate

A source can enter a canonical release only when:

1. checksums and provenance are complete;
2. sample and metadata alignment is explained;
3. genome build and coordinate convention are known;
4. critical warnings are resolved or explicitly accepted in an ADR;
5. aggregate reports are committed or archived;
6. conversion is deterministic and tested on a random subset.
