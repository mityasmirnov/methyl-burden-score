"""Orchestrate Stage 0 annotation graph builds."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from mbs.annotation.cgi_context import load_cpg_islands
from mbs.annotation.export_infinium import DEFAULT_PLATFORMS, load_infinium_probes
from mbs.annotation.gencode_regions import REGION_TYPES, ROLE_PRECEDENCE, build_gencode_regions
from mbs.annotation.locus_registry import build_locus_registry
from mbs.annotation.manifest import (
    git_commit,
    sha256_file,
    source_file_entry,
    utc_now_iso,
    validate_graph_manifest,
    write_json,
)
from mbs.annotation.map_loci import map_loci_to_regions, write_regions_bed
from mbs.annotation.regulatory_regions import build_rbs_regions
from mbs.annotation.tiles import DEFAULT_TILE_N, build_tiles
from mbs.paths import DataPaths

DEFAULT_GRAPH_ID = "graph-grch38-gencode38-five-role-v1"
GRAPH_V2_ID = "graph-grch38-gencode38-cgi-tile-v2"


@dataclass(frozen=True, slots=True)
class AnnotationGraphPaths:
    annotations_dir: Path
    graph_dir: Path
    report_dir: Path


def _default_report_dir(project_root: Path, graph_id: str) -> Path:
    slug = "annotation_graph_cgi_tile_v2" if graph_id == GRAPH_V2_ID else "annotation_graph_v1"
    return project_root / "reports" / "inspection" / slug


def default_paths(data_paths: DataPaths, graph_id: str = DEFAULT_GRAPH_ID) -> AnnotationGraphPaths:
    return AnnotationGraphPaths(
        annotations_dir=data_paths.data_root / "canonical" / "annotations",
        graph_dir=data_paths.data_root / "canonical" / "graphs" / graph_id,
        report_dir=_default_report_dir(data_paths.project_root, graph_id),
    )


def _write_parquet(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(path, index=False)


def _validation_report(
    *,
    loci: pd.DataFrame,
    probes: pd.DataFrame,
    edges: pd.DataFrame,
    genes: pd.DataFrame,
    regions: pd.DataFrame,
    locus_region_edges: pd.DataFrame,
    source_stats: dict[str, Any],
    region_policy: dict[str, Any],
) -> dict[str, Any]:
    if locus_region_edges.empty:
        assigned_loci: set[object] = set()
    else:
        assigned_loci = set(locus_region_edges["locus_id"].tolist())
    multi_gene = (
        locus_region_edges.merge(
            regions[["region_id", "gene_id"]],
            on="region_id",
            how="left",
        )
        .groupby("locus_id")["gene_id"]
        .nunique()
    )
    n_multi_gene = int((multi_gene > 1).sum()) if len(multi_gene) else 0
    probes_per_locus = pd.Series(dtype="int64") if edges.empty else edges.groupby("locus_id").size()
    region_counts = regions["region_type"].value_counts().to_dict() if not regions.empty else {}
    edge_role_counts = (
        locus_region_edges.merge(regions[["region_id", "region_type"]], on="region_id", how="left")[
            "region_type"
        ]
        .value_counts()
        .to_dict()
        if not locus_region_edges.empty
        else {}
    )
    platform_coverage = (
        probes.groupby("platform_id")
        .agg(
            n_probes=("probe_id", "size"),
            n_mapped=(
                "mapping_status",
                lambda s: int((s == "mapped").sum()),  # type: ignore[misc]
            ),
        )
        .reset_index()
        .to_dict(orient="records")
    )
    if "cpg_context" in loci.columns:
        cpg_context_counts = loci["cpg_context"].value_counts().astype(int).to_dict()
    else:
        cpg_context_counts = {}
    return {
        "source_stats": source_stats,
        "n_loci": len(loci),
        "n_genes": len(genes),
        "n_regions": len(regions),
        "n_locus_region_edges": len(locus_region_edges),
        "n_probes": len(probes),
        "n_probe_locus_edges": len(edges),
        "n_unmapped_probes": int((probes["mapping_status"] != "mapped").sum()),
        "n_unassigned_loci": int(len(loci) - len(assigned_loci)),
        "n_loci_multi_gene": n_multi_gene,
        "probes_collapsed_max": (
            int(probes_per_locus.max()) if len(probes_per_locus) else 0  # type: ignore[arg-type]
        ),
        "probes_collapsed_mean": (
            float(probes_per_locus.mean()) if len(probes_per_locus) else 0.0  # type: ignore[arg-type]
        ),
        "region_counts_by_type": {str(k): int(v) for k, v in region_counts.items()},
        "locus_region_edges_by_type": {str(k): int(v) for k, v in edge_role_counts.items()},
        "platform_coverage": platform_coverage,
        "cpg_context_counts": cpg_context_counts,
        "region_policy": region_policy,
    }


def attach_cgi_tile_systems(
    loci: pd.DataFrame,
    regions: pd.DataFrame,
    locus_region_edges: pd.DataFrame,
    *,
    tile_target_n_cpgs: int = DEFAULT_TILE_N,
    cgi_path: Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Append CGI RBS + CpG-count TBS; gene ∩ RBS empty; unmapped stay off tiles."""
    gene_regions = regions.copy()
    if "region_system" not in gene_regions.columns:
        gene_regions["region_system"] = "gene"
    assigned = set() if locus_region_edges.empty else set(locus_region_edges["locus_id"].tolist())
    islands = None
    if cgi_path is not None and cgi_path.is_file():
        islands = load_cpg_islands(cgi_path)
    rbs_r, rbs_e = build_rbs_regions(loci, gene_assigned_locus_ids=assigned, islands=islands)
    rbs_ids = set() if rbs_e.empty else set(rbs_e["locus_id"].tolist())
    leftover = set(loci["locus_id"].tolist()) - assigned - rbs_ids
    if "mapping_status" in loci.columns:
        mapped = set(loci.loc[loci["mapping_status"].eq("mapped"), "locus_id"].tolist())
        leftover &= mapped
    tbs_r, tbs_e = build_tiles(loci, remaining_locus_ids=leftover, target_n_cpgs=tile_target_n_cpgs)
    out_regions = pd.concat([gene_regions, rbs_r, tbs_r], ignore_index=True)
    out_edges = pd.concat([locus_region_edges, rbs_e, tbs_e], ignore_index=True)
    return out_regions, out_edges


def _v1_reuse_paths(data_root: Path) -> tuple[Path, Path, Path, Path, Path] | None:
    """Return v1 loci/probes/edges/genes/regions paths when all exist."""
    annotations_dir = data_root / "canonical" / "annotations"
    v1_graph = data_root / "canonical" / "graphs" / DEFAULT_GRAPH_ID
    loci_path = annotations_dir / "loci.parquet"
    probes_path = annotations_dir / "probes.parquet"
    probe_edges_path = annotations_dir / "probe_locus_edges.parquet"
    genes_path = v1_graph / "genes.parquet"
    regions_path = v1_graph / "regions.parquet"
    lr_path = v1_graph / "locus_region_edges.parquet"
    if not all(
        p.is_file()
        for p in (loci_path, probes_path, probe_edges_path, genes_path, regions_path, lr_path)
    ):
        return None
    return loci_path, probes_path, probe_edges_path, genes_path, regions_path


def build_annotation_graph(
    *,
    project_root: Path,
    data_root: Path,
    infinium_root: Path,
    gencode_path: Path,
    cgi_path: Path | None,
    graph_id: str = DEFAULT_GRAPH_ID,
    platforms: tuple[str, ...] | list[str] = DEFAULT_PLATFORMS,
    annotations_dir: Path | None = None,
    graph_dir: Path | None = None,
    report_dir: Path | None = None,
    probes: pd.DataFrame | None = None,
    genes: pd.DataFrame | None = None,
    regions: pd.DataFrame | None = None,
) -> dict[str, Any]:
    """Build locus registry + annotation graph artifacts under ``data_root``.

    For ``GRAPH_V2_ID``, reuse existing v1 annotations + five-role graph when
    present: write only the graph-v2 directory and v2 inspection report (do not
    clobber ``canonical/annotations/`` or ``annotation_graph_v1``).
    """
    annotations_dir = annotations_dir or (data_root / "canonical" / "annotations")
    graph_dir = graph_dir or (data_root / "canonical" / "graphs" / graph_id)
    report_dir = report_dir or _default_report_dir(project_root, graph_id)
    staging_dir = data_root / "staging" / "infinium_export"
    staging_dir.mkdir(parents=True, exist_ok=True)
    annotations_dir.mkdir(parents=True, exist_ok=True)
    graph_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)

    region_policy: dict[str, Any] = {
        "region_types": list(REGION_TYPES),
        "role_precedence": list(ROLE_PRECEDENCE),
        "coordinate_convention": "1-based_inclusive_cytosine; BED export 0-based_half_open",
        "intergenic_policy": "leave_unassigned",
        "promoter_core_bp": 200,
        "promoter_proximal_upstream_bp": 1500,
        "cgi_shore_bp": 2000,
        "cgi_shelf_bp": 4000,
        "gene_types": ["protein_coding"],
        "gencode_release": 38,
    }
    if graph_id == GRAPH_V2_ID:
        region_policy.update(
            {
                "region_systems": ["gene", "rbs", "tbs"],
                "tile_target_n_cpgs": DEFAULT_TILE_N,
                "rbs_source": "UCSC_cgi_per_island_shore",
            }
        )

    source_files: list[dict[str, Any]] = []
    source_stats: dict[str, Any] = {"platforms": list(platforms)}
    reuse_v1 = (
        graph_id == GRAPH_V2_ID
        and probes is None
        and genes is None
        and regions is None
        and _v1_reuse_paths(data_root) is not None
    )

    if reuse_v1:
        reuse = _v1_reuse_paths(data_root)
        if reuse is None:
            raise RuntimeError("v1 reuse paths vanished after check")
        loci_path, probes_path, probe_edges_path, v1_genes_path, v1_regions_path = reuse
        v1_lr_path = (
            data_root / "canonical" / "graphs" / DEFAULT_GRAPH_ID / "locus_region_edges.parquet"
        )
        loci = pd.read_parquet(loci_path)
        probes_out = pd.read_parquet(probes_path)
        probe_locus_edges = pd.read_parquet(probe_edges_path)
        genes = pd.read_parquet(v1_genes_path)
        regions = pd.read_parquet(v1_regions_path)
        locus_region_edges = pd.read_parquet(v1_lr_path)
        region_gene_edges = pd.read_parquet(
            data_root / "canonical" / "graphs" / DEFAULT_GRAPH_ID / "region_gene_edges.parquet"
        )
        source_stats["n_raw_probes"] = len(probes_out)
        source_stats["reused_v1_graph"] = DEFAULT_GRAPH_ID
        source_files.append(
            {
                "name": "reused_v1_annotations",
                "version": DEFAULT_GRAPH_ID,
                "uri": str(loci_path.resolve()),
                "sha256": sha256_file(loci_path),
                "license_note": "reuse existing Stage 0 annotations; not rewritten",
            }
        )
        source_files.append(
            {
                "name": "reused_v1_graph",
                "version": DEFAULT_GRAPH_ID,
                "uri": str(v1_genes_path.resolve()),
                "sha256": sha256_file(v1_genes_path),
                "license_note": "reuse five-role gene regions; attach RBS/TBS only",
            }
        )
        if cgi_path is not None and cgi_path.is_file():
            source_files.append(
                source_file_entry(
                    "ucsc_cpgIslandExt_hg38",
                    cgi_path,
                    version="UCSC_hg38",
                    uri=str(cgi_path),
                )
            )
        regions, locus_region_edges = attach_cgi_tile_systems(
            loci,
            regions,
            locus_region_edges,
            tile_target_n_cpgs=DEFAULT_TILE_N,
            cgi_path=cgi_path,
        )
        write_annotations = False
    else:
        if probes is None:
            probes_raw = load_infinium_probes(infinium_root, platforms)
            export_path = staging_dir / "probes_joined.parquet"
            _write_parquet(probes_raw, export_path)
            source_stats["n_raw_probes"] = len(probes_raw)
            for platform_id in platforms:
                ordering = infinium_root / platform_id / f"{platform_id}.ordering.tsv.gz"
                coord = infinium_root / platform_id / f"{platform_id}.hg38.coord.tsv.gz"
                source_files.append(
                    source_file_entry(
                        f"infinium_{platform_id}_ordering",
                        ordering,
                        version="InfiniumAnnotation_v8.1",
                        uri=str(ordering),
                    )
                )
                source_files.append(
                    source_file_entry(
                        f"infinium_{platform_id}_coord",
                        coord,
                        version="InfiniumAnnotation_v8.1",
                        uri=str(coord),
                    )
                )
        else:
            probes_raw = probes.copy()
            source_stats["n_raw_probes"] = len(probes_raw)
            source_files.append(
                {
                    "name": "fixture_probes",
                    "version": "test",
                    "uri": None,
                    "sha256": "0" * 64,
                    "license_note": "synthetic fixture",
                }
            )

        loci, probes_out, probe_locus_edges = build_locus_registry(probes_raw, cgi_path=cgi_path)

        if genes is None or regions is None:
            genes, regions = build_gencode_regions(gencode_path)
            source_files.append(
                source_file_entry(
                    "gencode_v38_annotation_gtf",
                    gencode_path,
                    version="GENCODE_v38",
                    uri=str(gencode_path),
                )
            )
        else:
            source_files.append(
                {
                    "name": "fixture_gencode",
                    "version": "test",
                    "uri": None,
                    "sha256": "0" * 64,
                    "license_note": "synthetic fixture",
                }
            )

        if cgi_path is not None and cgi_path.is_file():
            source_files.append(
                source_file_entry(
                    "ucsc_cpgIslandExt_hg38",
                    cgi_path,
                    version="UCSC_hg38",
                    uri=str(cgi_path),
                )
            )

        locus_region_edges, region_gene_edges = map_loci_to_regions(loci, regions)
        if graph_id == GRAPH_V2_ID:
            regions, locus_region_edges = attach_cgi_tile_systems(
                loci,
                regions,
                locus_region_edges,
                tile_target_n_cpgs=DEFAULT_TILE_N,
                cgi_path=cgi_path,
            )
        write_annotations = True
        loci_path = annotations_dir / "loci.parquet"
        probes_path = annotations_dir / "probes.parquet"
        probe_edges_path = annotations_dir / "probe_locus_edges.parquet"

    genes_path = graph_dir / "genes.parquet"
    regions_path = graph_dir / "regions.parquet"
    lr_path = graph_dir / "locus_region_edges.parquet"
    rg_path = graph_dir / "region_gene_edges.parquet"
    bed_path = graph_dir / "regions.bed"
    graph_manifest_path = graph_dir / "graph_manifest.json"
    annotations_manifest_path = annotations_dir / "annotations_manifest.json"
    validation_path = graph_dir / "validation_report.json"

    if write_annotations:
        _write_parquet(loci, loci_path)
        _write_parquet(probes_out, probes_path)
        _write_parquet(probe_locus_edges, probe_edges_path)
    _write_parquet(genes, genes_path)
    _write_parquet(regions, regions_path)
    _write_parquet(locus_region_edges, lr_path)
    _write_parquet(region_gene_edges, rg_path)
    write_regions_bed(regions, bed_path)

    builder_commit = git_commit(project_root)
    created_at = utc_now_iso()
    if graph_id == GRAPH_V2_ID:
        notes = (
            "Graph-v2: five-role gene + per-island CGI RBS + adaptive CpG-count TBS "
            "(ADR 0006). MethylCapsNet taxonomy reference only."
        )
    else:
        notes = (
            "Stage 0 five-role graph; MethylCapsNet used as taxonomy reference only; "
            "no CapsNet topology."
        )
    graph_manifest = {
        "artifact_version": "1",
        "graph_id": graph_id,
        "genome_build": "GRCh38",
        "builder_commit": builder_commit,
        "source_files": source_files,
        "region_policy": region_policy,
        "n_loci": len(loci),
        "n_genes": len(genes),
        "n_regions": len(regions),
        "n_locus_region_edges": len(locus_region_edges),
        "genes_path": str(genes_path.resolve()),
        "regions_path": str(regions_path.resolve()),
        "locus_region_edges_path": str(lr_path.resolve()),
        "region_gene_edges_path": str(rg_path.resolve()),
        "regions_bed_path": str(bed_path.resolve()),
        "created_at": created_at,
        "notes": notes,
    }
    validate_graph_manifest(graph_manifest)
    write_json(graph_manifest_path, graph_manifest)

    if write_annotations:
        annotations_manifest = {
            "artifact_version": "1",
            "genome_build": "GRCh38",
            "coordinate_convention": "1-based_cytosine_from_infinium_0based_CpG_beg",
            "builder_commit": builder_commit,
            "platforms": list(platforms),
            "source_files": source_files,
            "n_loci": len(loci),
            "n_probes": len(probes_out),
            "n_probe_locus_edges": len(probe_locus_edges),
            "loci_path": str(loci_path.resolve()),
            "probes_path": str(probes_path.resolve()),
            "probe_locus_edges_path": str(probe_edges_path.resolve()),
            "loci_sha256": sha256_file(loci_path),
            "created_at": created_at,
        }
        write_json(annotations_manifest_path, annotations_manifest)

    report = _validation_report(
        loci=loci,
        probes=probes_out,
        edges=probe_locus_edges,
        genes=genes,
        regions=regions,
        locus_region_edges=locus_region_edges,
        source_stats=source_stats,
        region_policy=region_policy,
    )
    write_json(validation_path, report)

    title = (
        "# Annotation graph validation (graph-v2 CGI + tiles)"
        if graph_id == GRAPH_V2_ID
        else "# Annotation graph validation (Stage 0)"
    )
    summary_md = report_dir / "summary.md"
    summary_md.write_text(
        "\n".join(
            [
                title,
                "",
                f"- graph_id: `{graph_id}`",
                f"- loci: {report['n_loci']}",
                f"- genes: {report['n_genes']}",
                f"- regions: {report['n_regions']}",
                f"- locus-region edges: {report['n_locus_region_edges']}",
                f"- unassigned/intergenic loci: {report['n_unassigned_loci']}",
                f"- multi-gene loci: {report['n_loci_multi_gene']}",
                f"- unmapped probes: {report['n_unmapped_probes']}",
                "",
                "## Region counts",
                "",
                *[
                    f"- `{name}`: {count}"
                    for name, count in sorted(report["region_counts_by_type"].items())
                ],
                "",
                "## Artifacts",
                "",
                f"- annotations: `{annotations_dir}`",
                f"- graph: `{graph_dir}`",
                f"- manifest: `{graph_manifest_path}`",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    write_json(report_dir / "summary.json", report)

    return {
        "graph_id": graph_id,
        "annotations_dir": str(annotations_dir),
        "graph_dir": str(graph_dir),
        "report_dir": str(report_dir),
        "graph_manifest_path": str(graph_manifest_path),
        "validation_report": report,
    }
