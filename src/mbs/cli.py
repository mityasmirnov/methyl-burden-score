"""Command-line interface for the MBS project."""

from __future__ import annotations

import json
import platform
import re
import shutil
import sys
from pathlib import Path
from typing import Annotated, Any

import torch
import typer
from rich.console import Console
from rich.table import Table

from mbs import __version__
from mbs.annotation.build import DEFAULT_GRAPH_ID, build_annotation_graph
from mbs.annotation.export_infinium import DEFAULT_PLATFORMS
from mbs.catalog import build_catalog, init_catalog
from mbs.inspect_cpgcorpus import inspect_cpgcorpus_gpl, write_cpgcorpus_report
from mbs.inspect_ewas_metadata import inspect_ewas_metadata, write_ewas_metadata_report
from mbs.inspect_source import inventory_source, write_inspection_report
from mbs.matrix.convert import DEFAULT_MATRIX_ID, convert_ewas_db_study
from mbs.paths import DataPaths, PathPolicyError
from mbs.static_features.export_cpgpt import DEFAULT_FEATURE_SET_ID, export_cpgpt_adapter
from mbs.training.loop import load_experiment_config, train_flat_baseline

app = typer.Typer(no_args_is_help=True, help="Methylation Burden Score tooling")
catalog_app = typer.Typer(no_args_is_help=True, help="DuckDB catalog operations")
inspect_app = typer.Typer(no_args_is_help=True, help="Source inspection reports")
graph_app = typer.Typer(no_args_is_help=True, help="Annotation graph builds")
matrix_app = typer.Typer(no_args_is_help=True, help="Canonical matrix conversion")
features_app = typer.Typer(no_args_is_help=True, help="Static locus feature export")
train_app = typer.Typer(no_args_is_help=True, help="Model training")
app.add_typer(catalog_app, name="catalog")
app.add_typer(inspect_app, name="inspect")
app.add_typer(graph_app, name="graph")
app.add_typer(matrix_app, name="matrix")
app.add_typer(features_app, name="features")
app.add_typer(train_app, name="train")
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


@inspect_app.command("ewas-metadata")
def inspect_ewas_metadata_cmd(
    report_dir: Annotated[
        Path | None,
        typer.Option(
            help="Report directory (default: reports/inspection/ewas_metadata_structure)"
        ),
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
            "logging": {"tensorboard": True},
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
        if study_holdout_fixture:
            cfg.setdefault("logging", {})["tensorboard"] = True

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
            }
        )
    )


if __name__ == "__main__":
    app()
