"""Command-line interface for the MBS project."""

from __future__ import annotations

import json
import platform
import re
import shutil
import sys
from pathlib import Path
from typing import Annotated, Any

import pandas as pd
import torch
import typer
from rich.console import Console
from rich.table import Table

from mbs import __version__
from mbs.annotation.build import DEFAULT_GRAPH_ID, GRAPH_V2_ID, build_annotation_graph
from mbs.annotation.export_infinium import DEFAULT_PLATFORMS
from mbs.catalog import build_catalog, init_catalog
from mbs.inspect_cpgcorpus import inspect_cpgcorpus_gpl, write_cpgcorpus_report
from mbs.inspect_ewas_metadata import inspect_ewas_metadata, write_ewas_metadata_report
from mbs.inspect_source import inventory_source, write_inspection_report
from mbs.matrix.convert import DEFAULT_MATRIX_ID, convert_ewas_db_study
from mbs.matrix.hub_pack import convert_hub_pack_subset, study_ids_from_sample_info
from mbs.matrix.hub_pack_index import (
    build_hub_pack_matrix_index,
    check_overlapping_gsm_betas,
)
from mbs.matrix.multitask_merge import merge_age_tissue_matrices
from mbs.matrix.store import read_sample_index
from mbs.matrix.virtual_hub_store import VIRTUAL_MATRIX_ID, build_virtual_hub_store
from mbs.paths import DataPaths, PathPolicyError
from mbs.release import (
    RELEASE_ID,
    refresh_release,
    release_paths,
    validate_release,
    write_phenotype_census_report,
    write_trait_eligibility_report,
)
from mbs.static_features.export_cpgpt import DEFAULT_FEATURE_SET_ID, export_cpgpt_adapter
from mbs.training.branch import train_branch_arm
from mbs.training.dev_cv import run_dev_cv
from mbs.training.hier_loop import train_hierarchical_baseline
from mbs.training.loop import load_experiment_config, train_flat_baseline
from mbs.training.monitor import run_monitor, ssh_tunnel_hint, validate_run_id
from mbs.training.phenotype_table import (
    HUB_UNION_PHENOTYPE_TABLE,
    HUB_UNION_SEX_ONTOLOGY,
    HUB_UNION_TISSUE_ONTOLOGY,
    build_hub_union_phenotype_table,
)

app = typer.Typer(no_args_is_help=True, help="Methylation Burden Score tooling")
catalog_app = typer.Typer(no_args_is_help=True, help="DuckDB catalog operations")
inspect_app = typer.Typer(no_args_is_help=True, help="Source inspection reports")
graph_app = typer.Typer(no_args_is_help=True, help="Annotation graph builds")
matrix_app = typer.Typer(no_args_is_help=True, help="Canonical matrix conversion")
features_app = typer.Typer(no_args_is_help=True, help="Static locus feature export")
train_app = typer.Typer(no_args_is_help=True, help="Model training")
phenotypes_app = typer.Typer(no_args_is_help=True, help="Phenotype table builds")
app.add_typer(catalog_app, name="catalog")
app.add_typer(inspect_app, name="inspect")
app.add_typer(graph_app, name="graph")
app.add_typer(matrix_app, name="matrix")
app.add_typer(features_app, name="features")
app.add_typer(train_app, name="train")
app.add_typer(phenotypes_app, name="phenotypes")
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


@catalog_app.command("refresh-release")
def catalog_refresh_release(
    release_id: Annotated[
        str,
        typer.Option(help="Versioned release id under canonical/releases/"),
    ] = RELEASE_ID,
    fetch_remote_index: Annotated[
        bool,
        typer.Option(
            "--fetch-remote-index/--no-fetch-remote-index",
            help="Optionally GET EWAS_db HTML index for advertised study count",
        ),
    ] = False,
    report_dir: Annotated[
        Path | None,
        typer.Option(
            help="Optional report dir (default: reports/inspection/{release_id})",
        ),
    ] = None,
) -> None:
    """Populate deepmat-data-v1 (or another release) from Hub Parquet + EWAS_db listing."""
    try:
        paths = DataPaths.from_environment()
    except PathPolicyError as error:
        console.print(f"[bold red]Path policy failure:[/bold red] {error}")
        raise typer.Exit(code=2) from error

    resolved_report = (
        report_dir.resolve()
        if report_dir is not None
        else paths.project_root / "reports" / "inspection" / release_id
    )
    _require_under_data(resolved_report, "report_dir")
    result = refresh_release(
        paths=paths,
        release_id=release_id,
        fetch_remote_index=fetch_remote_index,
        report_dir=resolved_report,
    )
    console.print_json(
        json.dumps(
            {
                "release_id": result.release_id,
                "release_root": result.release_root,
                "catalog_path": result.catalog_path,
                "n_samples": result.n_samples,
                "n_studies": result.n_studies,
                "n_phenotype_rows": result.n_phenotype_rows,
                "ewas_db_n_local_studies": result.ewas_db_n_local_studies,
                "ewas_db_n_local_gsm": result.ewas_db_n_local_gsm,
                "ewas_db_mirror_complete": result.ewas_db_mirror_complete,
                "report_dir": result.report_dir,
            }
        )
    )


@catalog_app.command("validate-release")
def catalog_validate_release(
    release_id: Annotated[str, typer.Option(help="Release id")] = RELEASE_ID,
) -> None:
    """Validate release_manifest.json and that the release catalog is populated."""
    try:
        paths = DataPaths.from_environment()
    except PathPolicyError as error:
        console.print(f"[bold red]Path policy failure:[/bold red] {error}")
        raise typer.Exit(code=2) from error
    try:
        result = validate_release(data_root=paths.data_root, release_id=release_id)
    except (FileNotFoundError, ValueError) as error:
        console.print(f"[bold red]validate-release failed:[/bold red] {error}")
        raise typer.Exit(code=1) from error
    console.print_json(json.dumps(result))


@catalog_app.command("phenotype-census")
def catalog_phenotype_census(
    release_id: Annotated[str, typer.Option(help="Release id")] = RELEASE_ID,
    report_dir: Annotated[
        Path | None,
        typer.Option(help="Report directory (default: reports/inspection/{release_id})"),
    ] = None,
) -> None:
    """Write phenotype census markdown/json from the release catalog."""
    try:
        paths = DataPaths.from_environment()
    except PathPolicyError as error:
        console.print(f"[bold red]Path policy failure:[/bold red] {error}")
        raise typer.Exit(code=2) from error
    rp = release_paths(paths.data_root, release_id)
    if not rp.catalog_db.is_file():
        console.print(f"[bold red]Release catalog missing:[/bold red] {rp.catalog_db}")
        raise typer.Exit(code=1)
    resolved_report = (
        report_dir.resolve()
        if report_dir is not None
        else paths.project_root / "reports" / "inspection" / release_id
    )
    _require_under_data(resolved_report, "report_dir")
    manifest = None
    if rp.manifest_path.is_file():
        manifest = json.loads(rp.manifest_path.read_text(encoding="utf-8"))
    out = write_phenotype_census_report(
        database=rp.catalog_db,
        report_dir=resolved_report,
        release_manifest=manifest,
    )
    console.print_json(json.dumps({"report_dir": str(out)}))


@catalog_app.command("trait-eligibility")
def catalog_trait_eligibility(
    release_id: Annotated[str, typer.Option(help="Release id")] = RELEASE_ID,
    report_dir: Annotated[
        Path | None,
        typer.Option(help="Report directory (default: reports/inspection/{release_id})"),
    ] = None,
) -> None:
    """Write trait eligibility report from the release catalog."""
    try:
        paths = DataPaths.from_environment()
    except PathPolicyError as error:
        console.print(f"[bold red]Path policy failure:[/bold red] {error}")
        raise typer.Exit(code=2) from error
    rp = release_paths(paths.data_root, release_id)
    if not rp.catalog_db.is_file():
        console.print(f"[bold red]Release catalog missing:[/bold red] {rp.catalog_db}")
        raise typer.Exit(code=1)
    resolved_report = (
        report_dir.resolve()
        if report_dir is not None
        else paths.project_root / "reports" / "inspection" / release_id
    )
    _require_under_data(resolved_report, "report_dir")
    out = write_trait_eligibility_report(database=rp.catalog_db, report_dir=resolved_report)
    console.print_json(json.dumps({"report_dir": str(out)}))


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


@inspect_app.command("ewas-metadata")
def inspect_ewas_metadata_cmd(
    report_dir: Annotated[
        Path | None,
        typer.Option(help="Report directory (default: reports/inspection/ewas_metadata_structure)"),
    ] = None,
) -> None:
    """Profile Atlas small tables and unpacked DataHub sample-info .txt packs."""
    try:
        paths = DataPaths.from_environment()
    except PathPolicyError as error:
        console.print(f"[bold red]Path policy failure:[/bold red] {error}")
        raise typer.Exit(code=2) from error

    resolved_report = report_dir or (
        paths.project_root / "reports" / "inspection" / "ewas_metadata_structure"
    )
    if not resolved_report.absolute().is_relative_to(paths.project_root.absolute()):
        resolved_report = _require_under_data(resolved_report, "report_dir")
    else:
        resolved_report = resolved_report.absolute()

    report = inspect_ewas_metadata(
        data_root=paths.data_root,
        project_root=paths.project_root,
    )
    written = write_ewas_metadata_report(report, resolved_report)
    present_packs = sum(1 for p in report["sample_packs"] if p.get("exists"))
    present_atlas = sum(1 for t in report["atlas_tables"] if t.get("exists"))
    console.print_json(
        json.dumps(
            {
                "report_dir": str(written),
                "atlas_tables_present": present_atlas,
                "sample_packs_present": present_packs,
                "generated_at": report["generated_at"],
            }
        )
    )


@graph_app.command("build")
def graph_build_cmd(
    graph_id: Annotated[
        str,
        typer.Option(help=f"Immutable graph release id (default five-role v1; also {GRAPH_V2_ID})"),
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


@matrix_app.command("convert")
def matrix_convert_cmd(  # noqa: PLR0917
    study_id: Annotated[
        str,
        typer.Option(help="Study accession (e.g. GSE35069)"),
    ] = "GSE35069",
    source_dir: Annotated[
        Path | None,
        typer.Option(help="EWAS_db study directory with GSM*.txt files"),
    ] = None,
    annotations_dir: Annotated[
        Path | None,
        typer.Option(help="Canonical annotations dir (loci + probe_locus_edges)"),
    ] = None,
    output_dir: Annotated[
        Path | None,
        typer.Option(help="Output matrix directory under canonical/matrices/"),
    ] = None,
    matrix_id: Annotated[
        str,
        typer.Option(help="Immutable matrix release id"),
    ] = DEFAULT_MATRIX_ID,
    platform_id: Annotated[
        str,
        typer.Option(help="Infinium platform id for probe→locus edges"),
    ] = "HM450",
    processing_level: Annotated[
        str,
        typer.Option(help="Processing level label (Hub baselines are GMQN)"),
    ] = "gmqn",
    report_dir: Annotated[
        Path | None,
        typer.Option(help="Inspection report directory (default under reports/inspection/)"),
    ] = None,
    verify: Annotated[
        bool,
        typer.Option(help="Round-trip verify raw Hub files against the written matrix"),
    ] = True,
) -> None:
    """Convert one EWAS Data Hub EWAS_db study into canonical matrix storage."""
    try:
        paths = DataPaths.from_environment()
    except PathPolicyError as error:
        console.print(f"[bold red]Path policy failure:[/bold red] {error}")
        raise typer.Exit(code=2) from error

    paths.ensure_directories()
    resolved_source = (
        source_dir or (paths.data_root / "raw" / "ewas_datahub" / "EWAS_db" / study_id)
    ).resolve()
    resolved_annotations = (
        annotations_dir or (paths.data_root / "canonical" / "annotations")
    ).resolve()
    resolved_output = (
        output_dir or (paths.data_root / "canonical" / "matrices" / matrix_id)
    ).resolve()
    resolved_report = (
        report_dir or (paths.project_root / "reports" / "inspection" / f"{study_id}_ewas_db")
    ).resolve()

    for label, path in (
        ("source_dir", resolved_source),
        ("annotations_dir", resolved_annotations),
        ("output_dir", resolved_output),
    ):
        _require_under_data(path, label)
    if not resolved_source.is_dir():
        raise typer.BadParameter(f"source_dir not found: {resolved_source}")
    if not (resolved_annotations / "probe_locus_edges.parquet").is_file():
        raise typer.BadParameter(
            f"annotations_dir missing probe_locus_edges.parquet: {resolved_annotations}"
        )

    result = convert_ewas_db_study(
        project_root=paths.project_root,
        source_dir=resolved_source,
        annotations_dir=resolved_annotations,
        output_dir=resolved_output,
        study_id=study_id,
        matrix_id=matrix_id,
        platform_id=platform_id,
        processing_level=processing_level,
        report_dir=resolved_report,
        verify=verify,
    )
    console.print_json(
        json.dumps(
            {
                "matrix_id": result.matrix_id,
                "study_id": result.study_id,
                "output_dir": str(result.output_dir),
                "report_dir": str(result.report_dir) if result.report_dir else None,
                "n_samples": result.stats["n_samples"],
                "n_study_loci": result.stats["n_study_loci"],
                "n_unmapped_probes": result.stats["n_unmapped_probes"],
                "roundtrip_ok": None if result.roundtrip is None else result.roundtrip.ok,
            }
        )
    )


@matrix_app.command("convert-pack")
def matrix_convert_pack_cmd(  # noqa: PLR0917
    phenotype_family: Annotated[
        str,
        typer.Option(
            help="Hub pack family: age|tissue|disease|cancer|blood|brain|sex|bmi|ancestry"
        ),
    ],
    matrix_id: Annotated[
        str,
        typer.Option(help="Immutable matrix release id"),
    ],
    study_ids: Annotated[
        str | None,
        typer.Option(help="Comma-separated study accessions (omit with --all-studies)"),
    ] = None,
    all_studies: Annotated[
        bool,
        typer.Option(
            "--all-studies",
            help="Include every study in the family's sample-info parquet",
        ),
    ] = False,
    annotations_dir: Annotated[
        Path | None,
        typer.Option(help="Canonical annotations dir (loci + probe_locus_edges)"),
    ] = None,
    output_dir: Annotated[
        Path | None,
        typer.Option(help="Output matrix directory under canonical/matrices/"),
    ] = None,
    sample_info: Annotated[
        Path | None,
        typer.Option(help="Sample-info parquet (default: canonical/phenotypes/<family>_…)"),
    ] = None,
    platform_id: Annotated[
        str,
        typer.Option(help="Infinium platform id for probe→locus edges"),
    ] = "HM450",
    processing_level: Annotated[
        str,
        typer.Option(help="Processing level label (Hub baselines are GMQN)"),
    ] = "gmqn",
    max_per_study: Annotated[
        int | None,
        typer.Option(help="Optional cap on samples per study (stable sample_id order)"),
    ] = None,
) -> None:
    """Convert a study-subset of one EWAS Data Hub baseline pack to a canonical matrix."""
    try:
        paths = DataPaths.from_environment()
    except PathPolicyError as error:
        console.print(f"[bold red]Path policy failure:[/bold red] {error}")
        raise typer.Exit(code=2) from error

    paths.ensure_directories()
    if all_studies and study_ids:
        raise typer.BadParameter("pass either --study-ids or --all-studies, not both")
    if not all_studies and not study_ids:
        raise typer.BadParameter("provide --study-ids or --all-studies")
    if max_per_study is not None and max_per_study < 1:
        raise typer.BadParameter("max_per_study must be >= 1")

    resolved_annotations = (
        annotations_dir or (paths.data_root / "canonical" / "annotations")
    ).resolve()
    resolved_output = (
        output_dir or (paths.data_root / "canonical" / "matrices" / matrix_id)
    ).resolve()
    for label, path in (
        ("annotations_dir", resolved_annotations),
        ("output_dir", resolved_output),
    ):
        _require_under_data(path, label)
    if not (resolved_annotations / "probe_locus_edges.parquet").is_file():
        raise typer.BadParameter(
            f"annotations_dir missing probe_locus_edges.parquet: {resolved_annotations}"
        )

    info_path = (
        sample_info.resolve()
        if sample_info is not None
        else (
            paths.data_root / "canonical" / "phenotypes" / f"{phenotype_family}_sample_info.parquet"
        )
    )
    if all_studies:
        if not info_path.is_file():
            raise typer.BadParameter(f"sample-info parquet required for --all-studies: {info_path}")
        studies = study_ids_from_sample_info(pd.read_parquet(info_path))
    else:
        if study_ids is None:
            raise typer.BadParameter("pass --study-ids or --all-studies")
        studies = [s.strip() for s in study_ids.split(",") if s.strip()]
        if not studies:
            raise typer.BadParameter("study_ids must list at least one accession")

    result = convert_hub_pack_subset(
        project_root=paths.project_root,
        data_root=paths.data_root,
        annotations_dir=resolved_annotations,
        phenotype_family=phenotype_family,
        study_ids=studies,
        matrix_id=matrix_id,
        output_dir=resolved_output,
        platform_id=platform_id,
        processing_level=processing_level,
        max_per_study=max_per_study,
        sample_info_path=info_path if sample_info is not None or all_studies else None,
    )
    console.print_json(
        json.dumps(
            {
                "matrix_id": result.matrix_id,
                "phenotype_family": result.phenotype_family,
                "study_ids": list(result.study_ids),
                "output_dir": str(result.output_dir),
                "n_samples": result.stats["n_samples"],
                "n_study_loci": result.stats["n_study_loci"],
                "n_unmapped_probes": result.stats["n_unmapped_probes"],
                "sample_phenotypes": result.stats["matrix_paths"]["sample_phenotypes"],
            }
        )
    )


@matrix_app.command("index-hub-packs")
def matrix_index_hub_packs_cmd(
    check_overlap: Annotated[
        bool,
        typer.Option(
            "--check-overlap/--no-check-overlap",
            help="Compare betas for GSMs in ≥2 pack matrices",
        ),
    ] = True,
    report_dir: Annotated[
        Path | None,
        typer.Option(help="Optional directory for overlap JSON report"),
    ] = None,
) -> None:
    """Build virtual Hub pack matrix index (+ optional GSM concordance check)."""
    try:
        paths = DataPaths.from_environment()
    except PathPolicyError as error:
        console.print(f"[bold red]Path policy failure:[/bold red] {error}")
        raise typer.Exit(code=2) from error

    paths.ensure_directories()
    index = build_hub_pack_matrix_index(paths.data_root)
    index_path = paths.data_root / "canonical" / "matrices" / "hub_pack_matrix_index.parquet"
    n_families = 0 if index.empty else len({str(x) for x in index["family"].tolist()})
    n_unique = 0 if index.empty else len({str(x) for x in index["sample_id"].tolist()})
    payload: dict[str, Any] = {
        "index_path": str(index_path),
        "n_rows": len(index),
        "n_families": n_families,
        "n_unique_gsm": n_unique,
    }
    if check_overlap:
        resolved_report = None
        if report_dir is not None:
            resolved_report_dir = _require_under_data(report_dir.resolve(), "report_dir")
            resolved_report = resolved_report_dir / "overlap_concordance.json"
            resolved_report.parent.mkdir(parents=True, exist_ok=True)
        overlap = check_overlapping_gsm_betas(
            paths.data_root,
            index=index,
            report_path=resolved_report,
        )
        payload["overlap"] = overlap.report
    console.print_json(json.dumps(payload))


@matrix_app.command("build-hub-virtual")
def matrix_build_hub_virtual_cmd(
    matrix_id: Annotated[
        str,
        typer.Option(help="Output virtual matrix id"),
    ] = VIRTUAL_MATRIX_ID,
) -> None:
    """Build Hub nine-pack virtual multi-store (route + indices; no dense Zarr)."""
    try:
        paths = DataPaths.from_environment()
    except PathPolicyError as error:
        console.print(f"[bold red]Path policy failure:[/bold red] {error}")
        raise typer.Exit(code=2) from error

    paths.ensure_directories()
    try:
        result = build_virtual_hub_store(
            data_root=paths.data_root,
            output_matrix_id=matrix_id,
        )
    except (FileNotFoundError, ValueError) as error:
        console.print(f"[bold red]Build failed:[/bold red] {error}")
        raise typer.Exit(code=1) from error

    console.print_json(
        json.dumps(
            {
                "matrix_id": result.matrix_id,
                "output_dir": str(result.output_dir),
                "n_samples": result.n_samples,
                "n_loci": result.n_loci,
                "n_source_matrices": result.n_source_matrices,
                "route_path": str(result.route_path),
                "stats": result.stats,
            },
            default=str,
        )
    )


@features_app.command("export-cpgpt")
def features_export_cpgpt_cmd(  # noqa: PLR0917
    feature_set_id: Annotated[
        str,
        typer.Option(help="Immutable static feature set id"),
    ] = DEFAULT_FEATURE_SET_ID,
    loci: Annotated[
        Path | None,
        typer.Option(help="Canonical loci.parquet (default: canonical/annotations/loci.parquet)"),
    ] = None,
    annotations_manifest: Annotated[
        Path | None,
        typer.Option(help="annotations_manifest.json for locus hash verification"),
    ] = None,
    output_dir: Annotated[
        Path | None,
        typer.Option(help="Output directory under canonical/static_features/"),
    ] = None,
    report_dir: Annotated[
        Path | None,
        typer.Option(help="Inspection report directory"),
    ] = None,
    device: Annotated[
        str,
        typer.Option(help="Torch device for encode_sequence (cuda or cpu)"),
    ] = "cuda",
    batch_size: Annotated[
        int,
        typer.Option(help="Adapter forward batch size"),
    ] = 8192,
) -> None:
    """Export CpGPT2M sequence-adapter embeddings for the locus registry."""
    try:
        paths = DataPaths.from_environment()
    except PathPolicyError as error:
        console.print(f"[bold red]Path policy failure:[/bold red] {error}")
        raise typer.Exit(code=2) from error

    paths.ensure_directories()
    annotations_dir = paths.data_root / "canonical" / "annotations"
    resolved_loci = (loci or (annotations_dir / "loci.parquet")).resolve()
    resolved_ann_manifest = (
        annotations_manifest or (annotations_dir / "annotations_manifest.json")
    ).resolve()
    resolved_output = (output_dir or paths.static_features_dir(feature_set_id)).resolve()
    resolved_report = (
        report_dir or (paths.project_root / "reports" / "inspection" / "static_features_cpgpt2m_v1")
    ).resolve()

    for label, path in (
        ("loci", resolved_loci.parent),
        ("output_dir", resolved_output),
    ):
        _require_under_data(path, label)
    if not resolved_loci.is_file():
        raise typer.BadParameter(f"loci parquet not found: {resolved_loci}")
    if not resolved_ann_manifest.is_file():
        raise typer.BadParameter(f"annotations_manifest.json not found: {resolved_ann_manifest}")
    if batch_size < 1:
        raise typer.BadParameter("batch_size must be >= 1")

    export_command = (
        f"uv run --extra cpgpt mbs features export-cpgpt --feature-set-id {feature_set_id}"
    )
    try:
        result = export_cpgpt_adapter(
            project_root=paths.project_root,
            loci_path=resolved_loci,
            annotations_manifest_path=resolved_ann_manifest,
            output_dir=resolved_output,
            report_dir=resolved_report,
            feature_set_id=feature_set_id,
            device=device,
            batch_size=batch_size,
            export_command=export_command,
        )
    except RuntimeError as error:
        console.print(f"[bold red]Export failed:[/bold red] {error}")
        raise typer.Exit(code=1) from error

    console.print_json(
        json.dumps(
            {
                "feature_set_id": result.feature_set_id,
                "output_dir": str(result.output_dir),
                "report_dir": str(result.report_dir) if result.report_dir else None,
                "n_loci": result.stats["n_loci"],
                "n_mapped": result.stats["n_mapped"],
                "n_missing": result.stats["n_missing"],
                "mapping_rate": result.stats["mapping_rate"],
                "checkpoint_sha256": result.manifest["checkpoint_sha256"],
                "locus_table_sha256": result.manifest["locus_table_sha256"],
            }
        )
    )


@phenotypes_app.command("build-multitask-table")
def phenotypes_build_multitask_table_cmd(  # noqa: PLR0917
    matrix_id: Annotated[
        str,
        typer.Option(help="Output multitask matrix id"),
    ] = "matrix-hub-age-tissue-multitask-v2",
    age_matrix_id: Annotated[
        str,
        typer.Option(help="Source age matrix id"),
    ] = "matrix-hub-age-studyholdout-v2",
    tissue_matrix_id: Annotated[
        str,
        typer.Option(help="Source tissue matrix id"),
    ] = "matrix-hub-tissue-studyholdout-v2",
    sex_matrix_id: Annotated[
        str | None,
        typer.Option(help="Optional source sex matrix id (Milestone 5d)"),
    ] = None,
    phenotype_table: Annotated[
        Path | None,
        typer.Option(help="Output sample phenotype parquet path"),
    ] = None,
    tissue_ontology: Annotated[
        Path | None,
        typer.Option(help="Output tissue ontology YAML path"),
    ] = None,
    sex_ontology: Annotated[
        Path | None,
        typer.Option(help="Output sex ontology YAML path"),
    ] = None,
    min_tissue_n: Annotated[
        int,
        typer.Option(help="Minimum samples per tissue class"),
    ] = 10,
) -> None:
    """Merge Hub matrices and write sample_phenotype_table + ontologies."""
    try:
        paths = DataPaths.from_environment()
    except PathPolicyError as error:
        console.print(f"[bold red]Path policy failure:[/bold red] {error}")
        raise typer.Exit(code=2) from error

    paths.ensure_directories()
    try:
        result = merge_age_tissue_matrices(
            project_root=paths.project_root,
            data_root=paths.data_root,
            age_matrix_id=age_matrix_id,
            tissue_matrix_id=tissue_matrix_id,
            sex_matrix_id=sex_matrix_id,
            output_matrix_id=matrix_id,
            phenotype_table_path=phenotype_table.resolve() if phenotype_table else None,
            tissue_ontology_path=tissue_ontology.resolve() if tissue_ontology else None,
            sex_ontology_path=sex_ontology.resolve() if sex_ontology else None,
            min_tissue_n=min_tissue_n,
        )
    except (FileNotFoundError, ValueError) as error:
        console.print(f"[bold red]Build failed:[/bold red] {error}")
        raise typer.Exit(code=1) from error

    console.print_json(
        json.dumps(
            {
                "matrix_id": result.matrix_id,
                "output_dir": str(result.output_dir),
                "n_samples": result.n_samples,
                "n_loci": result.n_loci,
                "n_deduped": result.n_deduped,
                "phenotype_table": str(result.phenotype_table_path),
                "tissue_ontology": str(result.tissue_ontology_path),
                "sex_ontology": (
                    None if result.sex_ontology_path is None else str(result.sex_ontology_path)
                ),
                "stats": result.stats,
            },
            default=str,
        )
    )


@phenotypes_app.command("build-hub-union-table")
def phenotypes_build_hub_union_table_cmd(
    matrix_id: Annotated[
        str,
        typer.Option(help="Virtual matrix id (sample/locus indices)"),
    ] = VIRTUAL_MATRIX_ID,
    phenotype_table: Annotated[
        Path | None,
        typer.Option(help="Output sample phenotype parquet path"),
    ] = None,
    tissue_ontology: Annotated[
        Path | None,
        typer.Option(help="Output tissue ontology YAML path"),
    ] = None,
    sex_ontology: Annotated[
        Path | None,
        typer.Option(help="Output sex ontology YAML path"),
    ] = None,
    min_tissue_n: Annotated[
        int,
        typer.Option(help="Minimum samples per tissue class"),
    ] = 10,
) -> None:
    """Join nine Hub pack phenotypes for the virtual multitask cohort (7E′)."""
    try:
        paths = DataPaths.from_environment()
    except PathPolicyError as error:
        console.print(f"[bold red]Path policy failure:[/bold red] {error}")
        raise typer.Exit(code=2) from error

    paths.ensure_directories()
    sample_index_path = (
        paths.data_root / "canonical" / "matrices" / matrix_id / "sample_index.parquet"
    )
    if not sample_index_path.is_file():
        console.print(
            f"[bold red]Missing virtual sample index:[/bold red] {sample_index_path}. "
            "Run `mbs matrix build-hub-virtual` first."
        )
        raise typer.Exit(code=1)
    sample_index = read_sample_index(sample_index_path)
    default_table = paths.data_root / "canonical" / "phenotypes" / HUB_UNION_PHENOTYPE_TABLE
    default_tissue = paths.data_root / "canonical" / "phenotypes" / HUB_UNION_TISSUE_ONTOLOGY
    default_sex = paths.data_root / "canonical" / "phenotypes" / HUB_UNION_SEX_ONTOLOGY
    try:
        result = build_hub_union_phenotype_table(
            data_root=paths.data_root,
            sample_index=sample_index,
            matrix_id=matrix_id,
            phenotype_table_path=phenotype_table.resolve() if phenotype_table else default_table,
            tissue_ontology_path=tissue_ontology.resolve() if tissue_ontology else default_tissue,
            sex_ontology_path=sex_ontology.resolve() if sex_ontology else default_sex,
            min_tissue_n=min_tissue_n,
        )
    except (FileNotFoundError, ValueError) as error:
        console.print(f"[bold red]Build failed:[/bold red] {error}")
        raise typer.Exit(code=1) from error

    console.print_json(
        json.dumps(
            {
                "matrix_id": matrix_id,
                "n_samples": result.n_samples,
                "phenotype_table": str(result.phenotype_table_path),
                "tissue_ontology": str(result.tissue_ontology_path),
                "sex_ontology": str(result.sex_ontology_path),
                "stats": result.stats,
            },
            default=str,
        )
    )


@train_app.command("flat")
def train_flat_cmd(  # noqa: PLR0917
    config: Annotated[
        Path | None,
        typer.Option(help="Experiment YAML (default: configs/experiment/stage0_flat_pilot.yaml)"),
    ] = None,
    run_id: Annotated[
        str,
        typer.Option(help="Artifact run id under $MBS_ARTIFACT_ROOT/runs/"),
    ] = "stage0-flat-gse35069-v1",
    device: Annotated[
        str,
        typer.Option(
            help="Torch device (cuda or cpu); cuda uses logical cuda:0 after CUDA_VISIBLE_DEVICES"
        ),
    ] = "cuda",
    overfit_fixture: Annotated[
        bool,
        typer.Option(
            "--overfit-fixture",
            help="Overfit the tiny synthetic fixture instead of the pilot matrix",
        ),
    ] = False,
    study_holdout_fixture: Annotated[
        bool,
        typer.Option(
            "--study-holdout-fixture",
            help="Multi-study synthetic fixture with study-grouped train/val/external_test",
        ),
    ] = False,
    max_epochs: Annotated[
        int | None,
        typer.Option(help="Override training.max_epochs"),
    ] = None,
    max_loci: Annotated[
        int | None,
        typer.Option(help="Optional study-column cap for smoke runs"),
    ] = None,
) -> None:
    """Train the exact flat DeepRVAT-style CpG→gene max-pooling baseline (deepMAT)."""
    try:
        paths = DataPaths.from_environment()
    except PathPolicyError as error:
        console.print(f"[bold red]Path policy failure:[/bold red] {error}")
        raise typer.Exit(code=2) from error

    paths.ensure_directories()
    default_config = paths.project_root / "configs" / "experiment" / "stage0_flat_pilot.yaml"
    config_path = (config or default_config).resolve()
    if overfit_fixture and study_holdout_fixture:
        raise typer.BadParameter("choose only one of --overfit-fixture / --study-holdout-fixture")
    fixture_mode = overfit_fixture or study_holdout_fixture
    if not fixture_mode and not config_path.is_file():
        raise typer.BadParameter(f"config not found: {config_path}")
    if max_epochs is not None and max_epochs < 1:
        raise typer.BadParameter("max_epochs must be >= 1")
    if max_loci is not None and max_loci < 1:
        raise typer.BadParameter("max_loci must be >= 1")
    if not run_id.strip() or "/" in run_id or ".." in run_id:
        raise typer.BadParameter("run_id must be a single path segment")

    if fixture_mode and not config_path.is_file():
        cfg: dict[str, Any] = {
            "experiment": {
                "name": (
                    "flat_study_holdout_fixture"
                    if study_holdout_fixture
                    else "flat_overfit_fixture"
                ),
                "stage": 0,
                "seed": 42,
            },
            "model": {
                "phi_layers": 2,
                "phi_hidden_dimension": 20,
                "rho_layers": 3,
                "rho_hidden_dimension": 10,
                "pooling": "max",
                "neutral_score": 0.5,
                "dropout": 0.0,
            },
            "training": {
                "optimizer": "adam",
                "learning_rate": 0.05,
                "weight_decay": 0.0,
                "max_epochs": max_epochs or 200,
                "early_stopping_patience": 50,
                "gradient_clip_norm": 2.0,
                "precision": "bf16-mixed",
                "require_cuda": False,
            },
            "logging": {"tensorboard": True, "auto_tensorboard": False},
            "heads": {"age": {"enabled": True}, "tissue": {"enabled": True}},
            "pilot": {"fixture_task": "tissue"},
        }
        if study_holdout_fixture:
            default_holdout = "stage0-flat-study-holdout-fixture"
            run_name = run_id if run_id != "stage0-flat-gse35069-v1" else default_holdout
        else:
            default_overfit = "stage0-flat-overfit-fixture"
            run_name = run_id if run_id != "stage0-flat-gse35069-v1" else default_overfit
    else:
        cfg = load_experiment_config(config_path)
        run_name = run_id
        if overfit_fixture or study_holdout_fixture:
            # Fixtures write TB events but skip spawning a server (CI / no training extra).
            log = cfg.setdefault("logging", {})
            log.setdefault("auto_tensorboard", False)
            if study_holdout_fixture:
                log["tensorboard"] = True

    try:
        result = train_flat_baseline(
            project_root=paths.project_root,
            data_root=paths.data_root,
            artifact_root=paths.artifact_root,
            config=cfg,
            run_id=run_name,
            device_str=device,
            overfit_fixture=overfit_fixture,
            study_holdout_fixture=study_holdout_fixture,
            max_epochs=max_epochs,
            max_loci=max_loci,
        )
    except RuntimeError as error:
        console.print(f"[bold red]Training failed:[/bold red] {error}")
        raise typer.Exit(code=1) from error

    console.print_json(
        json.dumps(
            {
                "run_id": result.run_id,
                "run_dir": str(result.run_dir),
                "checkpoint_dir": str(result.checkpoint_dir),
                "best_epoch": result.best_epoch,
                "final": result.metrics.get("final"),
                "overfit_ok": result.metrics.get("overfit_ok"),
                "device": result.metrics.get("device"),
                "tensorboard_url": result.tensorboard_url,
                "tensorboard_port": result.tensorboard_port,
                "monitor_hint": result.monitor_hint,
            }
        )
    )
    if result.tensorboard_url:
        console.print(
            f"[green]TensorBoard[/green] {result.tensorboard_url}  "
            f"({ssh_tunnel_hint(int(result.tensorboard_port or 6006))})"
        )
    if result.monitor_hint:
        console.print(f"[cyan]TUI monitor[/cyan] {result.monitor_hint}")


@train_app.command("hierarchical")
def train_hierarchical_cmd(  # noqa: PLR0917
    config: Annotated[
        Path | None,
        typer.Option(
            help="Experiment YAML (default: configs/experiment/stage0_hier_deeprvat_full.yaml)"
        ),
    ] = None,
    run_id: Annotated[
        str,
        typer.Option(help="Artifact run id under $MBS_ARTIFACT_ROOT/runs/"),
    ] = "stage0-hier-deeprvat-age-tissue-sex-full-v1",
    device: Annotated[
        str,
        typer.Option(
            help="Torch device (cuda or cpu); cuda uses logical cuda:0 after CUDA_VISIBLE_DEVICES"
        ),
    ] = "cuda",
    overfit_fixture: Annotated[
        bool,
        typer.Option(
            "--overfit-fixture",
            help="Overfit the tiny synthetic hierarchical fixture instead of the Hub matrix",
        ),
    ] = False,
    max_epochs: Annotated[
        int | None,
        typer.Option(help="Override training.max_epochs"),
    ] = None,
    max_loci: Annotated[
        int | None,
        typer.Option(help="Optional study-column cap for smoke runs"),
    ] = None,
) -> None:
    """Train HierarchicalDeepSet on the 5d age/tissue/sex cohort (Milestone 6)."""
    try:
        paths = DataPaths.from_environment()
    except PathPolicyError as error:
        console.print(f"[bold red]Path policy failure:[/bold red] {error}")
        raise typer.Exit(code=2) from error

    paths.ensure_directories()
    default_config = (
        paths.project_root / "configs" / "experiment" / "stage0_hier_deeprvat_full.yaml"
    )
    config_path = (config or default_config).resolve()
    if not overfit_fixture and not config_path.is_file():
        raise typer.BadParameter(f"config not found: {config_path}")
    if max_epochs is not None and max_epochs < 1:
        raise typer.BadParameter("max_epochs must be >= 1")
    if max_loci is not None and max_loci < 1:
        raise typer.BadParameter("max_loci must be >= 1")
    if not run_id.strip() or "/" in run_id or ".." in run_id:
        raise typer.BadParameter("run_id must be a single path segment")

    if overfit_fixture and not config_path.is_file():
        cfg: dict[str, Any] = {
            "experiment": {"name": "hier_overfit_fixture", "stage": 0, "seed": 42},
            "model": {
                "type": "hierarchical_deepset",
                "region_types": [
                    "promoter_core",
                    "promoter_proximal",
                    "five_prime",
                    "three_prime",
                    "gene_body",
                ],
                "residual_path": True,
                "cpg_hidden_dimension": 32,
                "region_hidden_dimension": 16,
                "region_type_dimension": 4,
                "cpg_pooling": "max",
                "region_pooling": "max",
                "residual_pooling": "max",
                "neutral_score": 0.5,
                "dropout": 0.0,
            },
            "training": {
                "optimizer": "adam",
                "learning_rate": 0.05,
                "weight_decay": 0.0,
                "max_epochs": max_epochs or 200,
                "early_stopping_patience": 50,
                "gradient_clip_norm": 2.0,
                "batch_size": 4,
                "precision": "bf16-mixed",
                "require_cuda": False,
            },
            "logging": {"tensorboard": False, "auto_tensorboard": False},
            "heads": {
                "age": {"enabled": True, "loss": "huber", "huber_delta": 1.0},
                "tissue": {"enabled": True},
                "sex": {"enabled": True},
            },
            "loss": {"lambda_age": 1.0, "lambda_tissue": 1.0, "lambda_sex": 1.0},
        }
    elif overfit_fixture:
        cfg = load_experiment_config(config_path)
        log = cfg.setdefault("logging", {})
        log.setdefault("auto_tensorboard", False)
    else:
        cfg = load_experiment_config(config_path)

    run_name = run_id
    if overfit_fixture and run_id == "stage0-hier-deeprvat-age-tissue-sex-full-v1":
        run_name = "stage0-hier-overfit-fixture"

    try:
        result = train_hierarchical_baseline(
            project_root=paths.project_root,
            data_root=paths.data_root,
            artifact_root=paths.artifact_root,
            config=cfg,
            run_id=run_name,
            device_str=device,
            overfit_fixture=overfit_fixture,
            max_epochs=max_epochs,
            max_loci=max_loci,
        )
    except RuntimeError as error:
        console.print(f"[bold red]Training failed:[/bold red] {error}")
        raise typer.Exit(code=1) from error

    console.print_json(
        json.dumps(
            {
                "run_id": result.run_id,
                "run_dir": str(result.run_dir),
                "checkpoint_dir": str(result.checkpoint_dir),
                "best_epoch": result.best_epoch,
                "final": result.metrics.get("final"),
                "external_test": result.metrics.get("external_test"),
                "ablations": result.metrics.get("ablations"),
                "annotation_slices": result.metrics.get("annotation_slices"),
                "n_residual_cols": result.metrics.get("n_residual_cols"),
                "vs_flat": result.metrics.get("vs_flat"),
                "overfit_ok": result.metrics.get("overfit_ok"),
                "device": result.metrics.get("device"),
                "tensorboard_url": result.tensorboard_url,
                "tensorboard_port": result.tensorboard_port,
                "monitor_hint": result.monitor_hint,
            },
            default=str,
        )
    )
    if result.tensorboard_url:
        console.print(
            f"[green]TensorBoard[/green] {result.tensorboard_url}  "
            f"({ssh_tunnel_hint(int(result.tensorboard_port or 6008))})"
        )
    if result.monitor_hint:
        console.print(f"[cyan]TUI monitor[/cyan] {result.monitor_hint}")


@train_app.command("branch")
def train_branch_cmd(  # noqa: PLR0917
    arm: Annotated[str, typer.Option(help="Independent arm: gene, rbs, tbs, or direct")] = "gene",
    run_id: Annotated[str, typer.Option(help="Artifact run id")] = "stage0-branch-fixture-v1",
    device: Annotated[str, typer.Option(help="Torch device")] = "cpu",
    config: Annotated[Path | None, typer.Option(help="Optional experiment YAML")] = None,
    overfit_fixture: Annotated[
        bool,
        typer.Option("--overfit-fixture/--hub", help="Synthetic fixture vs Hub config path"),
    ] = True,
    max_epochs: Annotated[int | None, typer.Option(help="Override training.max_epochs")] = None,
    max_loci: Annotated[int | None, typer.Option(help="Cap study loci for smoke")] = None,
) -> None:
    """Train one independently fitted branch arm (Milestone 7C/7E)."""
    try:
        paths = DataPaths.from_environment()
    except PathPolicyError as error:
        console.print(f"[bold red]Path policy failure:[/bold red] {error}")
        raise typer.Exit(code=2) from error
    paths.ensure_directories()
    if config is not None and config.is_file():
        cfg = load_experiment_config(config.resolve())
    else:
        cfg = {"experiment": {"seed": 42}, "model": {}, "training": {"max_epochs": 2}}

    result = train_branch_arm(
        arm=arm,
        project_root=paths.project_root,
        data_root=paths.data_root,
        artifact_root=paths.artifact_root,
        config=cfg,
        run_id=run_id,
        device=device,
        overfit_fixture=overfit_fixture,
        max_epochs=max_epochs,
        max_loci=max_loci,
    )
    console.print_json(json.dumps(result, default=str))


@train_app.command("cascade")
def train_cascade_cmd(  # noqa: PLR0917
    config: Annotated[
        Path,
        typer.Option(help="Experiment YAML (7F smoke or 7G eval)"),
    ] = Path("configs/experiment/stage0_7f_rbs_gene_direct.yaml"),
    run_id: Annotated[str, typer.Option(help="Artifact run id")] = "stage0-7f-cascade-v1",
    device: Annotated[str, typer.Option(help="Torch device")] = "cpu",
    overfit_fixture: Annotated[
        bool,
        typer.Option("--overfit-fixture/--hub", help="Synthetic fixture vs Hub frozen folds"),
    ] = False,
    max_folds: Annotated[
        int | None,
        typer.Option(help="Cap outer folds for Hub smoke (default: all)"),
    ] = None,
    max_loci: Annotated[
        int | None,
        typer.Option(help="Override cv_budget.max_loci for Hub smoke"),
    ] = None,
    max_epochs: Annotated[
        int | None,
        typer.Option(help="Override cv_budget.max_epochs for Hub smoke"),
    ] = None,
    max_train_samples: Annotated[
        int | None,
        typer.Option(help="Subsample train rows per fold (Hub smoke)"),
    ] = None,
    report_dir: Annotated[
        Path | None,
        typer.Option(help="Inspection report directory (default from YAML / milestone)"),
    ] = None,
    skip_if_done: Annotated[
        bool,
        typer.Option(
            "--skip-if-done/--force-retrain",
            help="Skip folds that already have scores/score_manifest.json",
        ),
    ] = True,
) -> None:
    """RBS→gene cascade + leftover direct (no TBS). 7F smoke or 7G Hub eval."""
    from mbs.training.cascade_loop import train_cascade
    from mbs.training.loop import load_experiment_config

    try:
        paths = DataPaths.from_environment()
    except PathPolicyError as error:
        console.print(f"[bold red]Path policy failure:[/bold red] {error}")
        raise typer.Exit(code=2) from error
    paths.ensure_directories()
    cfg_path = config if config.is_absolute() else paths.project_root / config
    cfg = None if overfit_fixture else load_experiment_config(cfg_path)
    if cfg is not None:
        budget = dict(cfg.get("cv_budget") or {})
        if max_loci is not None:
            budget["max_loci"] = int(max_loci)
        if max_epochs is not None:
            budget["max_epochs"] = int(max_epochs)
        cfg["cv_budget"] = budget
    resolved_report: Path | None = None
    if report_dir is not None:
        resolved_report = report_dir if report_dir.is_absolute() else paths.project_root / report_dir
        resolved_report = _require_under_data(resolved_report.resolve(), "report_dir")
    result = train_cascade(
        project_root=paths.project_root,
        data_root=paths.data_root,
        artifact_root=paths.artifact_root,
        config=cfg,
        config_path=None if overfit_fixture else cfg_path,
        run_id=run_id,
        device_str=device,
        overfit_fixture=overfit_fixture,
        max_folds=max_folds,
        max_train_samples=max_train_samples,
        report_dir=resolved_report,
        skip_if_done=skip_if_done,
    )
    console.print_json(json.dumps(result.metrics, default=str))
    console.print(f"[green]report[/green] {result.report_dir}")


@train_app.command("dev-cv")
def train_dev_cv_cmd(  # noqa: PLR0917
    config: Annotated[
        Path,
        typer.Option(help="Bake-off YAML (default: stage0_7e_bakeoff.yaml)"),
    ] = Path("configs/experiment/stage0_7e_bakeoff.yaml"),
    report_dir: Annotated[
        Path,
        typer.Option(help="Inspection report directory"),
    ] = Path("reports/inspection/stage0_7e_dev_cv"),
    device: Annotated[str, typer.Option(help="Torch device")] = "cuda",
    max_epochs: Annotated[int | None, typer.Option(help="Override max_epochs for all arms")] = None,
    max_loci: Annotated[int | None, typer.Option(help="Cap loci (Hub smoke / CV budget)")] = None,
    fixture: Annotated[
        bool,
        typer.Option("--fixture/--hub", help="Synthetic fold plumbing only (not Done when)"),
    ] = False,
) -> None:
    """Milestone 7E development CV orchestrator (3×2 arms, shared folds)."""
    try:
        paths = DataPaths.from_environment()
    except PathPolicyError as error:
        console.print(f"[bold red]Path policy failure:[/bold red] {error}")
        raise typer.Exit(code=2) from error
    paths.ensure_directories()
    cfg_path = config if config.is_absolute() else paths.project_root / config
    out_dir = report_dir if report_dir.is_absolute() else paths.project_root / report_dir
    summary = run_dev_cv(
        project_root=paths.project_root,
        data_root=paths.data_root,
        artifact_root=paths.artifact_root,
        bakeoff_config=cfg_path,
        report_dir=out_dir,
        device=device,
        max_epochs=max_epochs,
        max_loci=max_loci,
        fixture=fixture,
    )
    console.print_json(json.dumps(summary, default=str))


@app.command("monitor")
def monitor_cmd(  # noqa: PLR0917
    run_id: Annotated[
        str,
        typer.Option(help="Run id under $MBS_ARTIFACT_ROOT/runs/<run_id>/"),
    ],
    config: Annotated[
        Path | None,
        typer.Option(help="Experiment YAML for max_epochs when resolved_config is absent"),
    ] = None,
    max_epochs: Annotated[
        int | None,
        typer.Option(help="Override training.max_epochs for ETA"),
    ] = None,
    interval: Annotated[
        float,
        typer.Option(help="Refresh interval in seconds"),
    ] = 2.0,
    once: Annotated[
        bool,
        typer.Option("--once", help="Print one snapshot and exit (no live refresh)"),
    ] = False,
    no_tensorboard: Annotated[
        bool,
        typer.Option(
            "--no-tensorboard",
            help="Do not start/reuse TensorBoard (TUI only)",
        ),
    ] = False,
    tb_port: Annotated[
        int,
        typer.Option(help="Preferred TensorBoard port (falls back if busy)"),
    ] = 6006,
) -> None:
    """Live TUI + TensorBoard for a train run (default: both)."""
    try:
        paths = DataPaths.from_environment()
    except PathPolicyError as error:
        console.print(f"[bold red]Path policy failure:[/bold red] {error}")
        raise typer.Exit(code=2) from error

    try:
        validate_run_id(run_id)
    except ValueError as error:
        raise typer.BadParameter(str(error)) from error

    config_path: Path | None = None
    if config is not None:
        config_path = _require_under_data(config.absolute(), "config")
        if not config_path.is_file():
            raise typer.BadParameter(f"config not found: {config_path}")
    if max_epochs is not None and max_epochs < 1:
        raise typer.BadParameter("max_epochs must be >= 1")
    if interval < 0.2:
        raise typer.BadParameter("interval must be >= 0.2")
    if tb_port < 1 or tb_port > 65535:
        raise typer.BadParameter("tb_port out of range")

    try:
        snap = run_monitor(
            run_id=run_id,
            artifact_root=paths.artifact_root,
            config_path=config_path,
            max_epochs=max_epochs,
            interval_s=interval,
            once=once,
            start_tensorboard=not no_tensorboard,
            tb_port=tb_port,
        )
    except ValueError as error:
        console.print(f"[bold red]Monitor error:[/bold red] {error}")
        raise typer.Exit(code=1) from error
    except RuntimeError as error:
        console.print(f"[bold red]Monitor error:[/bold red] {error}")
        raise typer.Exit(code=1) from error
    except KeyboardInterrupt:
        console.print("\n[dim]monitor stopped[/dim]")
        raise typer.Exit(code=0) from None

    if snap.tensorboard is not None:
        console.print(
            f"[green]TensorBoard[/green] {snap.tensorboard.url}  "
            f"({'reused' if snap.tensorboard.reused else 'started'})  "
            f"tunnel: {ssh_tunnel_hint(snap.tensorboard.port)}"
        )
    if snap.status == "finished":
        console.print(f"[green]Run finished:[/green] {snap.run_root}")
    elif snap.status == "stalled":
        raise typer.Exit(code=1)


if __name__ == "__main__":
    app()
