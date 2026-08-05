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

Keep DataHub external validation releases separate from CpGCorpus training releases until study overlap and normalization provenance are resolved.

Inspect:

- study and GSE identifiers;
- trait bundle definitions;
- sample information file;
- matrix processing and batch-correction description;
- platform and genome build;
- overlap with CpGCorpus studies;
- phenotype coding.

## EWAS Atlas

Atlas is an association evidence source, not a sample matrix. Freeze and checksum association, study, cohort, and annotation exports. Map CpGs through the project graph rather than relying only on supplied nearest-gene fields.

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
