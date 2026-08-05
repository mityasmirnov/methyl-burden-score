"""Command-line interface for the MBS project."""

from __future__ import annotations

import json
import platform
import shutil
import sys
from pathlib import Path
from typing import Annotated

import torch
import typer
from rich.console import Console
from rich.table import Table

from mbs import __version__
from mbs.catalog import build_catalog
from mbs.paths import DataPaths, PathPolicyError

app = typer.Typer(no_args_is_help=True, help="Methylation Burden Score tooling")
catalog_app = typer.Typer(no_args_is_help=True, help="DuckDB catalog operations")
app.add_typer(catalog_app, name="catalog")
console = Console()


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
        if not path.absolute().is_relative_to(Path("/data")):
            raise typer.BadParameter(f"path must be under /data: {path}")

    result = build_catalog(
        database=database,
        sql_dir=sql_dir,
        parquet_root=parquet_root,
        read_only=read_only,
    )
    console.print_json(json.dumps(result))


if __name__ == "__main__":
    app()
