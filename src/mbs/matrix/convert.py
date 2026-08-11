"""Convert one EWAS Data Hub study into canonical matrix storage."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from mbs.annotation.manifest import git_commit, sha256_file, utc_now_iso, write_json
from mbs.matrix.ewas_db import (
    beta_qc_stats,
    list_ewas_db_sample_files,
    read_ewas_db_sample,
)
from mbs.matrix.locus_map import build_probe_locus_map, load_probe_locus_edges
from mbs.matrix.roundtrip import RoundTripResult, verify_roundtrip
from mbs.matrix.store import (
    ARTIFACT_VERSION,
    DEFAULT_DTYPE,
    MATRIX_ORIENTATION,
    MISSING_VALUE_ENCODING,
    matrix_store_paths,
    source_file_records,
    write_betas_zarr,
    write_locus_index,
    write_matrix_manifest,
    write_sample_index,
)

DEFAULT_MATRIX_ID = "matrix-gse35069-ewasdb-v1"


@dataclass(frozen=True, slots=True)
class ConvertResult:
    matrix_id: str
    study_id: str
    output_dir: Path
    report_dir: Path | None
    stats: dict[str, Any]
    roundtrip: RoundTripResult | None


def convert_ewas_db_study(
    *,
    project_root: Path,
    source_dir: Path,
    annotations_dir: Path,
    output_dir: Path,
    study_id: str,
    matrix_id: str,
    platform_id: str = "HM450",
    processing_level: str = "gmqn",
    report_dir: Path | None = None,
    verify: bool = True,
    verify_max_probes: int = 2048,
    verify_sample_ids: list[str] | None = None,
) -> ConvertResult:
    """Write one Hub study as ``betas.zarr`` + indices + ``matrix_manifest.json``."""
    source_dir = source_dir.resolve()
    annotations_dir = annotations_dir.resolve()
    output_dir = output_dir.resolve()
    project_root = project_root.resolve()

    sample_files = list_ewas_db_sample_files(source_dir)
    sample_ids = [s.sample_id for s in sample_files]
    source_sample_ids = [s.source_sample_id for s in sample_files]

    # Read all samples once; build sorted union probe vocabulary.
    tables = [
        read_ewas_db_sample(sample.path, sample_id=sample.sample_id) for sample in sample_files
    ]
    probe_set: set[str] = set()
    for table in tables:
        probe_set.update(str(p) for p in table.probe_ids)
    observed_probes = np.asarray(sorted(probe_set), dtype=object)

    edges = load_probe_locus_edges(annotations_dir, platform_id=platform_id)
    locus_map = build_probe_locus_map(observed_probes, edges, platform_id=platform_id)

    n_samples = len(sample_files)
    n_loci = len(locus_map.locus_ids)
    betas = np.full((n_samples, n_loci), np.nan, dtype=np.float32)
    probe_to_col = {str(pid): i for i, pid in enumerate(locus_map.probe_ids)}

    aggregate_finite: list[np.ndarray] = []
    n_out_of_range = 0
    n_missing_cells = 0

    for row, table in enumerate(tables):
        mapped_cols = pd.Series(table.probe_ids.astype(str)).map(probe_to_col)
        valid = mapped_cols.notna().to_numpy()
        if valid.any():
            cols = mapped_cols.to_numpy(dtype=np.float64)[valid].astype(np.int64)
            values = table.betas.astype(np.float32, copy=False)[valid]
            betas[row, cols] = values
        row_qc = beta_qc_stats(betas[row, :])
        n_out_of_range += int(row_qc["n_out_of_range"])
        n_missing_cells += int(row_qc["n_missing"])
        finite_mask = np.isfinite(betas[row, :])
        if finite_mask.any():
            aggregate_finite.append(betas[row, finite_mask].copy())

    paths = matrix_store_paths(output_dir)
    paths.root.mkdir(parents=True, exist_ok=True)
    write_sample_index(
        paths.sample_index_path,
        sample_ids=sample_ids,
        source_sample_ids=source_sample_ids,
    )
    write_locus_index(
        paths.locus_index_path,
        locus_ids=locus_map.locus_ids,
        canonical_keys=locus_map.canonical_keys,
        probe_ids=locus_map.probe_ids,
        annotation_status=locus_map.annotation_status,
    )
    chunks = (min(64, max(1, n_samples)), min(4096, max(1, n_loci)))
    write_betas_zarr(paths.betas_path, betas, chunks=chunks)

    source_records = source_file_records([s.path for s in sample_files])
    sample_index_sha = sha256_file(paths.sample_index_path)
    locus_index_sha = sha256_file(paths.locus_index_path)

    if aggregate_finite:
        all_finite = np.concatenate(aggregate_finite)
        beta_min = float(all_finite.min())
        beta_max = float(all_finite.max())
        beta_mean = float(all_finite.mean())
        n_finite = int(all_finite.size)
    else:
        beta_min = float("nan")
        beta_max = float("nan")
        beta_mean = float("nan")
        n_finite = 0

    notes = (
        f"EWAS Data Hub EWAS_db/{study_id}; platform={platform_id}; "
        f"unmapped_probes={len(locus_map.unmapped_probe_ids)}; "
        f"residual_probes={locus_map.n_residual_probes}; "
        f"collapsed_probes={locus_map.n_collapsed_probes}"
    )
    manifest: dict[str, Any] = {
        "artifact_version": ARTIFACT_VERSION,
        "matrix_id": matrix_id,
        "study_id": study_id,
        "platform_id": platform_id,
        "processing_level": processing_level,
        "genome_build": "GRCh38",
        "shape": [n_samples, n_loci],
        "dtype": DEFAULT_DTYPE,
        "chunks": list(chunks),
        "compression": None,
        "missing_value_encoding": MISSING_VALUE_ENCODING,
        "matrix_orientation": MATRIX_ORIENTATION,
        "matrix_path": str(paths.betas_path),
        "sample_index_path": str(paths.sample_index_path),
        "locus_index_path": str(paths.locus_index_path),
        "sample_index_sha256": sample_index_sha,
        "locus_index_sha256": locus_index_sha,
        "source_files": source_records,
        "conversion_commit": git_commit(project_root),
        "created_at": utc_now_iso(),
        "notes": notes,
    }
    write_matrix_manifest(paths.manifest_path, manifest)

    stats: dict[str, Any] = {
        "n_samples": n_samples,
        "n_study_loci": n_loci,
        "n_observed_probes": locus_map.n_observed_probes,
        "n_mapped_probes": locus_map.n_mapped_probes,
        "n_unmapped_probes": len(locus_map.unmapped_probe_ids),
        "n_residual_probes": locus_map.n_residual_probes,
        "n_collapsed_probes": locus_map.n_collapsed_probes,
        "n_finite_betas": n_finite,
        "n_missing_cells": n_missing_cells,
        "n_out_of_range": n_out_of_range,
        "beta_min": beta_min,
        "beta_max": beta_max,
        "beta_mean": beta_mean,
        "unmapped_probe_ids_head": list(locus_map.unmapped_probe_ids[:50]),
        "matrix_paths": {
            "root": str(paths.root),
            "betas": str(paths.betas_path),
            "sample_index": str(paths.sample_index_path),
            "locus_index": str(paths.locus_index_path),
            "manifest": str(paths.manifest_path),
        },
    }

    roundtrip: RoundTripResult | None = None
    if verify:
        check_samples = verify_sample_ids
        if check_samples is None:
            # First, middle, last sample when available
            idxs = sorted({0, n_samples // 2, n_samples - 1})
            check_samples = [sample_ids[i] for i in idxs]
        roundtrip = verify_roundtrip(
            source_dir,
            output_dir,
            sample_ids=check_samples,
            max_probes=verify_max_probes,
        )
        stats["roundtrip"] = {
            "ok": roundtrip.ok,
            "n_compared": roundtrip.n_compared,
            "n_mismatch": roundtrip.n_mismatch,
            "max_abs_diff": roundtrip.max_abs_diff,
            "sample_ids": list(roundtrip.sample_ids),
        }
        if not roundtrip.ok:
            raise RuntimeError(
                f"round-trip verification failed: mismatches={roundtrip.n_mismatch} "
                f"max_abs_diff={roundtrip.max_abs_diff}"
            )

    written_report: Path | None = None
    if report_dir is not None:
        written_report = write_conversion_report(report_dir, stats=stats, manifest=manifest)

    return ConvertResult(
        matrix_id=matrix_id,
        study_id=study_id,
        output_dir=paths.root,
        report_dir=written_report,
        stats=stats,
        roundtrip=roundtrip,
    )


def write_conversion_report(
    report_dir: Path,
    *,
    stats: dict[str, Any],
    manifest: dict[str, Any],
) -> Path:
    report_dir = report_dir.resolve()
    report_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "study_id": manifest["study_id"],
        "matrix_id": manifest["matrix_id"],
        "platform_id": manifest["platform_id"],
        "processing_level": manifest["processing_level"],
        "genome_build": manifest["genome_build"],
        "shape": manifest["shape"],
        "dtype": manifest["dtype"],
        "conversion_commit": manifest["conversion_commit"],
        "created_at": manifest["created_at"],
        "stats": stats,
        "source_file_count": len(manifest["source_files"]),
        "source_files_head": manifest["source_files"][:5],
    }
    write_json(report_dir / "summary.json", payload)
    lines = [
        f"# Inspection: {manifest['study_id']} EWAS Data Hub → canonical matrix",
        "",
        f"- matrix_id: `{manifest['matrix_id']}`",
        f"- platform_id: `{manifest['platform_id']}`",
        f"- processing_level: `{manifest['processing_level']}`",
        f"- genome_build: `{manifest['genome_build']}`",
        f"- shape: `{manifest['shape']}` (samples x loci)",
        f"- dtype: `{manifest['dtype']}`",
        f"- conversion_commit: `{manifest['conversion_commit']}`",
        "",
        "## Stats",
        "",
        "```json",
        json.dumps(stats, indent=2, sort_keys=True),
        "```",
        "",
        "## Round-trip",
        "",
    ]
    rt = stats.get("roundtrip")
    if rt is None:
        lines.append("Not run.")
    elif rt.get("ok"):
        lines.append(
            f"PASS — compared {rt['n_compared']} cells across "
            f"{len(rt['sample_ids'])} samples; max_abs_diff={rt['max_abs_diff']}."
        )
    else:
        lines.append(f"FAIL — mismatches={rt['n_mismatch']} max_abs_diff={rt['max_abs_diff']}.")
    lines.append("")
    (report_dir / "summary.md").write_text("\n".join(lines), encoding="utf-8")
    return report_dir
