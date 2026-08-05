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
from mbs.catalog import build_catalog, init_catalog
from mbs.inspect_source import inventory_source, write_inspection_report
from mbs.paths import DataPaths, PathPolicyError

app = typer.Typer(no_args_is_help=True, help="Methylation Burden Score tooling")
catalog_app = typer.Typer(no_args_is_help=True, help="DuckDB catalog operations")
inspect_app = typer.Typer(no_args_is_help=True, help="Source inspection reports")
app.add_typer(catalog_app, name="catalog")
app.add_typer(inspect_app, name="inspect")
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


if __name__ == "__main__":
    app()
