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

Milestone **7C** graph-v2 layers ([ADR 0006](adr/0006-multipath-noncoding-scores.md)):

```text
locus -> non-gene regulatory region (cCRE / enhancer / CGI / DMR / ChromHMM)  → RBS
locus -> adaptive CpG-count intergenic tile                                    → TBS
remaining / explicit residual loci                                             → direct CpG path
```

Do **not** use one catch-all residual capsule for all leftover CpGs (MethylCapsNet
`include_last`-style). Prefer multiple biologically or spatially meaningful
region systems; shared Deep Sets scale better than one free network per capsule.

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

Stage 0 five-role graph (`graph-grch38-gencode38-five-role-v1`) records them as
unassigned canonical loci. Hierarchical-v0.1 routes them to a one-scalar
residual path (frozen baseline only).

**Graph v2 / Milestone 7C on-disk assignment order** (tile nodes may remain):

1. Gene-linked five-role regions;
2. Non-gene regulatory annotations (cCRE, enhancer, CGI/shore, DMR, ChromHMM,
   DHS) without requiring a gene → **RBS**;
3. Remaining coordinate-mapped loci → adaptive **CpG-count tiles** (graph
   nodes only after 7F);
4. Isolated / unmapped / explicitly retained loci → **direct CpG** contribution.

**Score topology from Milestone 7F** ([ADR 0009](adr/0009-drop-tbs-scores.md)):
typed regions (gene roles + RBS) → region scores; allocate RBS to genes
(typed role and/or nearest-gene of the **region**) → **MBS**; orphan RBS stay
genome-wide; leftover CpGs (including former tile-only loci) → **direct**.
Do not emit TBS scores. Do not nearest-gene leftover **CpGs**
([ADR 0004](adr/0004-unmapped-probe-retention.md)).

Any regulatory or tile edge requires an evidence type, source version, and
confidence field. Optional later: cCRE-to-gene or eQTM-supported gene edges
with the same provenance fields—never silent nearest-gene of unmapped CpGs.

### MethylCapsNet lessons

Use: multiple capsule / region systems (genes, CGI, enhancers, bins, custom
BED). Avoid: dropping CpGs absent from a selected capsule; dropping small
capsules solely for size; one free network per capsule at genome scale; one
global “all remaining CpGs” module as the production path.

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

Changing an interval boundary, source annotation, overlap precedence, or edge
policy requires a new graph release. Graph-v2 (RBS/TBS) uses the immutable
identifier `graph-grch38-gencode38-cgi-tile-v2` (per-island CGI RBS + adaptive
CpG-count TBS; ADR 0006). Build via
`mbs graph build --graph-id graph-grch38-gencode38-cgi-tile-v2` (reuses v1
annotations; writes only the graph-v2 directory). **On disk** under
`$MBS_DATA_ROOT/canonical/graphs/graph-grch38-gencode38-cgi-tile-v2/`
(inspection: `reports/inspection/annotation_graph_cgi_tile_v2/`; RBS
`UCSC_cgi_per_island_shore`). Plan:
[`plans/milestone-7c-graph-v2-topology.md`](plans/milestone-7c-graph-v2-topology.md).