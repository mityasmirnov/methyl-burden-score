# Annotation coverage (probe → locus → five-role region)

- graph_id: `graph-grch38-gencode38-five-role-v1`
- loci: 1082522
- assigned loci: 761410 (70.34%)
- unassigned loci: 321112 (29.66%)
- multi-gene loci: 92793
- locus-region edges: 870947
- probes (all platforms): 2290670
- unmapped probes: 2143 (0.09%)

## Figures

![Locus assigned vs unassigned](figures/locus_assigned_pie.png)

![Per-array mapped probe assignment](figures/platform_assigned_vs_unassigned.png)

![Loci by regulatory role](figures/loci_by_role.png)

![CpG-island context](figures/cpg_island_context.png)

## Definitions

- **Unmapped probe:** Illumina probe with no GRCh38 cytosine coordinate; excluded from matrices
- **Unassigned locus:** Mapped locus with no five-role GENCODE region edge; hier train uses singleton unassigned
- **Atlas:** EWAS Atlas probe TSV is a separate unused layer; not counted here

## Locus-level roles (unique loci with ≥1 edge of that type)

- `promoter_core`: 187602 loci (200605 edges)
- `promoter_proximal`: 169085 loci (175070 edges)
- `five_prime`: 21967 loci (22071 edges)
- `gene_body`: 403607 loci (433506 edges)
- `three_prime`: 38526 loci (39695 edges)

## Per-platform probe coverage

| Platform | Probes | Mapped | Unmapped | Mapped→region | Mapped unassigned |
|----------|-------:|-------:|---------:|--------------:|------------------:|
| EPIC | 866553 | 865904 | 649 (0.0749%) | 624518 (72.1232%) | 241386 (27.8768%) |
| EPICv2 | 937690 | 937054 | 636 (0.0678%) | 668862 (71.3792%) | 268192 (28.6208%) |
| HM450 | 486427 | 485569 | 858 (0.1764%) | 375826 (77.3991%) | 109743 (22.6009%) |

## Per-platform probes by role (unique probes)

### `EPIC`

- `promoter_core`: 167912
- `promoter_proximal`: 146326
- `five_prime`: 18668
- `gene_body`: 312246
- `three_prime`: 31902

### `EPICv2`

- `promoter_core`: 164300
- `promoter_proximal`: 148370
- `five_prime`: 18730
- `gene_body`: 355831
- `three_prime`: 33452

### `HM450`

- `promoter_core`: 115163
- `promoter_proximal`: 91481
- `five_prime`: 15186
- `gene_body`: 166646
- `three_prime`: 24293

## Cross-check vs `annotation_graph_v1`

- ok: `True`
- unassigned loci here / ref: 321112 / 321112

## Regenerate

```bash
uv run python scripts/write_annotation_coverage_report.py
uv sync --extra analysis  # once, for matplotlib
uv run python scripts/write_pipeline_doc_figures.py
```
