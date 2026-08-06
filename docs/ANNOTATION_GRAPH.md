# Annotation graph

## Purpose

The annotation graph separates biological topology from model code. It maps canonical GRCh38 CpG loci to typed regulatory regions and regions to genes.

Stage 0 uses only evidence that can be built reproducibly and audited across 450K, EPIC, and EPICv2 arrays.

Implementation plan and rebuild commands for Milestone 2:
[`plans/milestone-2-canonical-annotation-graph.md`](plans/milestone-2-canonical-annotation-graph.md).
CLI: `mbs graph build`.

## Graph layers

```text
probe -> canonical locus -> gene region -> gene
```

Future-compatible optional layers:

```text
locus -> cCRE or enhancer capsule
locus -> fixed genomic tile
locus -> pathway capsule
```

## Canonical sources

Initial sources should be pinned by checksum and version:

- GENCODE release 38 protein-coding genes and transcripts;
- Illumina 450K, EPIC, and EPICv2 manifests remapped to GRCh38;
  - reannotation code: `vendor/epicv2_manifest` ([EPICv2_manifest](https://github.com/bethan-mallabar-rimmer/EPICv2_manifest));
  - published tables: Zenodo ([doi:10.5281/zenodo.14933468](https://doi.org/10.5281/zenodo.14933468)), stored under `$MBS_DATA_ROOT/raw/manifests/`;
- UCSC CpG-island intervals;
- a frozen consensus cCRE or ChromHMM annotation release;
- probe-quality and cross-reactivity annotations where licensing permits redistribution.

MethylCapsNet resources are used to define useful annotation categories, not copied as the canonical graph. Its legacy implementation is primarily hg19 and 450K-oriented.

## Gene-region roles

Stage 0 creates five mutually prioritized gene roles:

1. `promoter_core`: TSS -200 to TSS +200;
2. `promoter_proximal`: TSS -1500 to TSS -200;
3. `five_prime`: first exon and 5-prime UTR not assigned above;
4. `three_prime`: 3-prime UTR not assigned above;
5. `gene_body`: transcribed span not assigned to a higher-priority role.

The interval definitions are configurable and recorded in the graph manifest.

## Transcript consolidation

Transcript annotations are retained in source tables, but Stage 0 builds gene-level union regions.

For one locus–gene pair, primary-role precedence is:

```text
promoter_core
> promoter_proximal
> five_prime
> three_prime
> gene_body
```

This avoids counting the same CpG repeatedly because of overlapping transcripts.

A locus may still map to multiple genes when gene intervals genuinely overlap. Those edges remain explicit.

## Orthogonal locus annotations

Each locus or locus–region edge may carry:

### CpG context

```text
island
north_shore
south_shore
north_shelf
south_shelf
open_sea
```

### Regulatory context

```text
promoter-like cCRE
enhancer-like cCRE
CTCF-only cCRE
DHS-only cCRE
ChromHMM state
DHS flag
```

### Probe context

```text
array platforms
probe design
cross-reactivity flag
SNP-overlap flag
quality flag
```

These are features or QC fields; they do not create an uncontrolled cross-product of region types.

## Intergenic loci

Do not force intergenic loci to the nearest gene.

Stage 0 records them as unassigned canonical loci. Later graph releases may add:

- cCRE-to-gene evidence edges;
- eQTM-supported edges;
- fixed genomic tiles;
- user-defined capsules.

Any such edge requires an evidence type, source version, and confidence field.

## Graph tables

### `genes.parquet`

```text
gene_id
gene_name
chromosome
start
end
strand
gene_type
source_version
```

### `regions.parquet`

```text
region_id
gene_id
region_type
chromosome
start
end
strand
source_version
```

### `locus_region_edges.parquet`

```text
locus_id
region_id
edge_weight
evidence_type
primary_gene_role
```

### `region_gene_edges.parquet`

```text
region_id
gene_id
edge_weight
```

## Validation report

Every graph build reports:

- input file hashes and row counts;
- loci, regions, genes, and edges;
- region counts by type;
- loci per region and per gene;
- loci assigned to more than one gene;
- probes collapsed to one locus;
- ambiguous and unmapped probes;
- unassigned/intergenic loci;
- coverage by platform, study, gene, and region type;
- coordinate and contig mismatches;
- interval-policy parameters.

## BED export

The graph builder exports the same regulatory regions to BED:

```text
chrom start end region_id score strand gene_id region_type
```

This BED is the future interface to methylartist and other long-read aggregation tools, ensuring that array and long-read validation use identical region boundaries.

## Versioning

A graph release is immutable after use in an experiment. A semantic graph identifier should encode the genome build and policy version, for example:

```text
graph-grch38-gencode38-five-role-v1
```

Changing an interval boundary, source annotation, overlap precedence, or edge policy requires a new graph release.