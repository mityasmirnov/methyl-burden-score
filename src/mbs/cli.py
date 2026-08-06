"""Command-line interface for the MBS project."""

from __future__ import annotations

import json
import platform
import re
import shutil
import sys
from pathlib import Path
from typing import Annotated

import torch
import typer
from rich.console import Console
from rich.table import Table

from mbs import __version__
from mbs.annotation.build import DEFAULT_GRAPH_ID, build_annotation_graph
from mbs.annotation.export_infinium import DEFAULT_PLATFORMS
from mbs.catalog import build_catalog, init_catalog
from mbs.inspect_cpgcorpus import inspect_cpgcorpus_gpl, write_cpgcorpus_report
from mbs.inspect_source import inventory_source, write_inspection_report
from mbs.paths import DataPaths, PathPolicyError

app = typer.Typer(no_args_is_help=True, help="Methylation Burden Score tooling")
catalog_app = typer.Typer(no_args_is_help=True, help="DuckDB catalog operations")
inspect_app = typer.Typer(no_args_is_help=True, help="Source inspection reports")
graph_app = typer.Typer(no_args_is_help=True, help="Annotation graph builds")
app.add_typer(catalog_app, name="catalog")
app.add_typer(inspect_app, name="inspect")
app.add_typer(graph_app, name="graph")
console = Console()

_SOURCE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")


def _require_under_data(path: Path, label: str) -> Path:
    absolute = path.absolute()
    if not absolute.is_relative_to(Path("/data")):
        raise typer.BadParameter(f"{label} must be under /data: {absolute}")
    return absolute


@app.command()
def doctor(create_directories: bool = False) -> None:
    """Validate the data-only workspace and print environment diagnostics."""
    try:
        paths = DataPaths.from_environment()
    except PathPolicyError as error:
        console.print(f"[bold red]Path policy failure:[/bold red] {error}")
        raise typer.Exit(code=2) from error

    if create_directories:
        paths.ensure_directories()

    table = Table(title="MBS environment")
    table.add_column("Item")
    table.add_column("Value")
    table.add_column("Status")

    table.add_row("MBS version", __version__, "ok")
    table.add_row("Python", sys.version.split()[0], "ok")
    table.add_row("Platform", platform.platform(), "ok")
    table.add_row("PyTorch", torch.__version__, "ok")
    table.add_row("CUDA available", str(torch.cuda.is_available()), "ok")
    if torch.cuda.is_available():
        table.add_row("CUDA devices", str(torch.cuda.device_count()), "ok")

    for name, raw_path in paths.as_dict().items():
        path = Path(raw_path)
        exists = path.exists()
        status = "ok" if exists else "missing"
        table.add_row(name, raw_path, status)

    usage = shutil.disk_usage("/data")
    free_gib = usage.free / (1024**3)
    table.add_row("/data free", f"{free_gib:,.1f} GiB", "ok" if free_gib > 20 else "low")

    console.print(table)

    if not paths.project_root.exists():
        console.print("[yellow]Project root is missing; run scripts/bootstrap_server.sh.[/yellow]")


@catalog_app.command("init")
def catalog_init(
    database: Annotated[
        Path | None,
        typer.Option(
            help="Output DuckDB path (default: $MBS_DATA_ROOT/canonical/catalog/catalog.duckdb)"
        ),
    ] = None,
    sql_dir: Annotated[
        Path | None,
        typer.Option(help="Directory containing numbered SQL files (default: $MBS_ROOT/sql)"),
    ] = None,
    parquet_root: Annotated[
        Path | None,
        typer.Option(help="Parquet table directory (default: .../canonical/catalog/tables)"),
    ] = None,
) -> None:
    """Create Stage 0 directories and apply the catalog SQL schema."""
    try:
        paths = DataPaths.from_environment()
    except PathPolicyError as error:
        console.print(f"[bold red]Path policy failure:[/bold red] {error}")
        raise typer.Exit(code=2) from error

    for label, path in (
        ("database", database),
        ("sql_dir", sql_dir),
        ("parquet_root", parquet_root),
    ):
        if path is not None:
            _require_under_data(path, label)

    result = init_catalog(
        paths=paths,
        sql_dir=sql_dir,
        database=database,
        parquet_root=parquet_root,
    )
    console.print_json(json.dumps(result))


@catalog_app.command("build")
def catalog_build(
    database: Annotated[Path, typer.Option(help="Output DuckDB path")],
    sql_dir: Annotated[Path, typer.Option(help="Directory containing numbered SQL files")],
    parquet_root: Annotated[
        Path,
        typer.Option(help="Canonical Parquet table directory exposed to SQL as parquet_root"),
    ],
    read_only: bool = False,
) -> None:
    """Build or validate the analytical catalog."""
    paths = DataPaths.from_environment()
    paths.validate()

    for path in (database, sql_dir, parquet_root):
        _require_under_data(path, "path")

    result = build_catalog(
        database=database,
        sql_dir=sql_dir,
        parquet_root=parquet_root,
        read_only=read_only,
    )
    console.print_json(json.dumps(result))


@inspect_app.command("source")
def inspect_source_cmd(
    source_id: Annotated[str, typer.Option(help="Source identifier, e.g. cpgcorpus")],
    raw_root: Annotated[
        Path | None,
        typer.Option(help="Raw source directory (default: $MBS_DATA_ROOT/raw/{source_id})"),
    ] = None,
    report_dir: Annotated[
        Path | None,
        typer.Option(help="Report directory (default: reports/inspection/{source_id})"),
    ] = None,
    max_entries: Annotated[int, typer.Option(help="Maximum file entries to inventory")] = 10_000,
) -> None:
    """Write a shallow inventory report for one raw source directory."""
    if not _SOURCE_ID_RE.fullmatch(source_id):
        raise typer.BadParameter(
            "source_id must be 1-128 chars: letters, digits, underscore, dot, or hyphen"
        )
    if max_entries < 1:
        raise typer.BadParameter("max_entries must be >= 1")

    try:
        paths = DataPaths.from_environment()
    except PathPolicyError as error:
        console.print(f"[bold red]Path policy failure:[/bold red] {error}")
        raise typer.Exit(code=2) from error

    resolved_raw = _require_under_data(
        raw_root or (paths.data_root / "raw" / source_id),
        "raw_root",
    )
    resolved_report = report_dir or (paths.project_root / "reports" / "inspection" / source_id)
    if not resolved_report.absolute().is_relative_to(paths.project_root.absolute()):
        resolved_report = _require_under_data(resolved_report, "report_dir")
    else:
        resolved_report = resolved_report.absolute()

    inventory = inventory_source(
        resolved_raw,
        source_id=source_id,
        max_entries=max_entries,
    )
    written = write_inspection_report(inventory, resolved_report)
    console.print_json(
        json.dumps(
            {
                "source_id": source_id,
                "raw_root": str(resolved_raw),
                "report_dir": str(written),
                "file_count": inventory["file_count"],
                "total_bytes": inventory["total_bytes"],
                "truncated": inventory["truncated"],
            }
        )
    )


@inspect_app.command("cpgcorpus-gpl")
def inspect_cpgcorpus_gpl_cmd(
    gse: Annotated[str, typer.Option(help="GSE accession, e.g. GSE125367")],
    gpl: Annotated[str, typer.Option(help="GPL platform, e.g. GPL21145")],
    raw_root: Annotated[
        Path | None,
        typer.Option(help="GSE/GPL directory (default: raw/cpgcorpus/{gse}/{gpl})"),
    ] = None,
    report_dir: Annotated[
        Path | None,
        typer.Option(help="Report directory (default: reports/inspection/{gse}_{gpl})"),
    ] = None,
) -> None:
    """Inspect one CpGCorpus GSE/GPL for layout, alignment, betas, and metadata."""
    if not _SOURCE_ID_RE.fullmatch(gse) or not _SOURCE_ID_RE.fullmatch(gpl):
        raise typer.BadParameter("gse and gpl must be safe accession-like identifiers")

    try:
        paths = DataPaths.from_environment()
    except PathPolicyError as error:
        console.print(f"[bold red]Path policy failure:[/bold red] {error}")
        raise typer.Exit(code=2) from error

    resolved_raw = _require_under_data(
        raw_root or (paths.data_root / "raw" / "cpgcorpus" / gse / gpl),
        "raw_root",
    )
    source_id = f"{gse}_{gpl}"
    resolved_report = report_dir or (paths.project_root / "reports" / "inspection" / source_id)
    if not resolved_report.absolute().is_relative_to(paths.project_root.absolute()):
        resolved_report = _require_under_data(resolved_report, "report_dir")
    else:
        resolved_report = resolved_report.absolute()

    report = inspect_cpgcorpus_gpl(resolved_raw, gse=gse, gpl=gpl)
    written = write_cpgcorpus_report(report, resolved_report)
    console.print_json(
        json.dumps(
            {
                "source_id": report["source_id"],
                "report_dir": str(written),
                "perfect_alignment": report["sample_alignment"].get("perfect_alignment"),
                "n_samples": report["value_qc"].get("n_samples"),
                "n_probes": report["value_qc"].get("n_probes"),
                "warnings": report["warnings"],
            }
        )
    )


@graph_app.command("build")
def graph_build_cmd(
    graph_id: Annotated[
        str,
        typer.Option(help="Immutable graph release id"),
    ] = DEFAULT_GRAPH_ID,
    platforms: Annotated[
        str,
        typer.Option(help="Comma-separated InfiniumAnnotation platform dirs"),
    ] = ",".join(DEFAULT_PLATFORMS),
    infinium_root: Annotated[
        Path | None,
        typer.Option(
            help="InfiniumAnnotation root (default: $MBS_ROOT/vendor/infinium_annotation)"
        ),
    ] = None,
    gencode: Annotated[
        Path | None,
        typer.Option(help="GENCODE GTF path (default: raw/gencode/gencode.v38.annotation.gtf.gz)"),
    ] = None,
    cgi: Annotated[
        Path | None,
        typer.Option(
            help="UCSC CpG island table (default: raw/annotations/cpgIslandExt.hg38.txt.gz)"
        ),
    ] = None,
) -> None:
    """Build the Stage 0 locus registry and five-role annotation graph."""
    try:
        paths = DataPaths.from_environment()
    except PathPolicyError as error:
        console.print(f"[bold red]Path policy failure:[/bold red] {error}")
        raise typer.Exit(code=2) from error

    paths.ensure_directories()
    default_infinium = paths.project_root / "vendor" / "infinium_annotation"
    resolved_infinium = (infinium_root or default_infinium).resolve()
    resolved_gencode = (
        gencode or (paths.data_root / "raw" / "gencode" / "gencode.v38.annotation.gtf.gz")
    ).resolve()
    resolved_cgi = (
        cgi or (paths.data_root / "raw" / "annotations" / "cpgIslandExt.hg38.txt.gz")
    ).resolve()

    _require_under_data(paths.data_root, "data_root")
    if not resolved_gencode.is_file():
        raise typer.BadParameter(f"GENCODE GTF not found: {resolved_gencode}")
    if not resolved_infinium.is_dir():
        raise typer.BadParameter(f"InfiniumAnnotation root not found: {resolved_infinium}")
    cgi_path = resolved_cgi if resolved_cgi.is_file() else None

    platform_list = tuple(p.strip() for p in platforms.split(",") if p.strip())
    if not platform_list:
        raise typer.BadParameter("platforms must list at least one platform id")

    result = build_annotation_graph(
        project_root=paths.project_root,
        data_root=paths.data_root,
        infinium_root=resolved_infinium,
        gencode_path=resolved_gencode,
        cgi_path=cgi_path,
        graph_id=graph_id,
        platforms=platform_list,
    )
    console.print_json(
        json.dumps(
            {
                "graph_id": result["graph_id"],
                "annotations_dir": result["annotations_dir"],
                "graph_dir": result["graph_dir"],
                "report_dir": result["report_dir"],
                "n_loci": result["validation_report"]["n_loci"],
                "n_genes": result["validation_report"]["n_genes"],
                "n_regions": result["validation_report"]["n_regions"],
                "n_locus_region_edges": result["validation_report"]["n_locus_region_edges"],
            }
        )
    )


if __name__ == "__main__":
    app()
