"""Canonical matrix store (Zarr betas + parquet indices + manifest)."""

from __future__ import annotations

import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import zarr

from mbs.annotation.manifest import sha256_file, write_json

_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_COMMIT_RE = re.compile(r"^[0-9a-f]{40}$")
_DATA_PATH_RE = re.compile(r"^/data/")

ARTIFACT_VERSION = "matrix-manifest-v1"
MISSING_VALUE_ENCODING = "nan"
MATRIX_ORIENTATION = "samples_by_loci"
DEFAULT_DTYPE = "float32"


@dataclass(frozen=True, slots=True)
class MatrixStorePaths:
    root: Path
    betas_path: Path
    sample_index_path: Path
    locus_index_path: Path
    manifest_path: Path


def matrix_store_paths(output_dir: Path) -> MatrixStorePaths:
    root = output_dir.resolve()
    return MatrixStorePaths(
        root=root,
        betas_path=root / "betas.zarr",
        sample_index_path=root / "sample_index.parquet",
        locus_index_path=root / "locus_index.parquet",
        manifest_path=root / "matrix_manifest.json",
    )


def validate_matrix_manifest(manifest: dict[str, Any]) -> None:
    """Validate required fields against ``schemas/matrix_manifest.schema.json`` rules."""
    required = [
        "artifact_version",
        "matrix_id",
        "study_id",
        "platform_id",
        "processing_level",
        "genome_build",
        "shape",
        "dtype",
        "missing_value_encoding",
        "matrix_path",
        "sample_index_path",
        "locus_index_path",
        "source_files",
        "conversion_commit",
        "created_at",
    ]
    missing = [key for key in required if key not in manifest]
    if missing:
        raise ValueError(f"matrix manifest missing keys: {missing}")
    if manifest["genome_build"] != "GRCh38":
        raise ValueError("genome_build must be GRCh38")
    if manifest.get("matrix_orientation", MATRIX_ORIENTATION) != MATRIX_ORIENTATION:
        raise ValueError("matrix_orientation must be samples_by_loci")
    if manifest["dtype"] not in {"float16", "float32", "float64"}:
        raise ValueError(f"unsupported dtype: {manifest['dtype']}")
    shape = manifest["shape"]
    if not isinstance(shape, list) or len(shape) != 2:
        raise ValueError("shape must be [n_samples, n_loci]")
    if any(not isinstance(x, int) or x < 0 for x in shape):
        raise ValueError("shape entries must be non-negative integers")
    if not _COMMIT_RE.fullmatch(str(manifest["conversion_commit"])):
        raise ValueError("conversion_commit must be a 40-char lowercase hex SHA")
    for key in ("matrix_path", "sample_index_path", "locus_index_path"):
        path = str(manifest[key])
        if not _DATA_PATH_RE.match(path):
            raise ValueError(f"{key} must be an absolute /data path, got {path}")
    source_files = manifest["source_files"]
    if not isinstance(source_files, list) or not source_files:
        raise ValueError("source_files must be a non-empty list")
    for entry in source_files:
        if not isinstance(entry, dict):
            raise TypeError("source_files entries must be objects")
        if "path" not in entry or "sha256" not in entry:
            raise ValueError("source_files entry requires path and sha256")
        if not _DATA_PATH_RE.match(str(entry["path"])):
            raise ValueError(f"source path must be under /data/: {entry['path']}")
        if not _SHA256_RE.fullmatch(str(entry["sha256"])):
            raise ValueError(f"invalid sha256 for {entry['path']}")
    for key in ("sample_index_sha256", "locus_index_sha256"):
        value = manifest.get(key)
        if value is not None and not _SHA256_RE.fullmatch(str(value)):
            raise ValueError(f"invalid {key}")


def write_sample_index(
    path: Path,
    *,
    sample_ids: list[str],
    source_sample_ids: list[str],
) -> pd.DataFrame:
    if len(sample_ids) != len(source_sample_ids):
        raise ValueError("sample_ids and source_sample_ids length mismatch")
    frame = pd.DataFrame(
        {
            "row_index": np.arange(len(sample_ids), dtype=np.int64),
            "sample_id": sample_ids,
            "source_sample_id": source_sample_ids,
        }
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(path, index=False)
    return frame


def write_locus_index(
    path: Path,
    *,
    locus_ids: np.ndarray,
    canonical_keys: np.ndarray,
    probe_ids: np.ndarray,
) -> pd.DataFrame:
    n = len(locus_ids)
    if len(canonical_keys) != n or len(probe_ids) != n:
        raise ValueError("locus index arrays must have equal length")
    frame = pd.DataFrame(
        {
            "col_index": np.arange(n, dtype=np.int64),
            "locus_id": locus_ids.astype(np.uint64, copy=False),
            "canonical_key": canonical_keys.astype(str),
            "probe_id": probe_ids.astype(str),
        }
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(path, index=False)
    return frame


def write_betas_zarr(
    path: Path,
    betas: np.ndarray,
    *,
    chunks: tuple[int, int] | None = None,
) -> tuple[int, int]:
    """Write ``[n_samples, n_loci]`` float32 betas to Zarr (NaN = missing)."""
    if betas.ndim != 2:
        raise ValueError(f"betas must be 2-D, got shape {betas.shape}")
    if betas.dtype != np.float32:
        betas = betas.astype(np.float32, copy=False)
    if path.exists():
        shutil.rmtree(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    n_samples, n_loci = int(betas.shape[0]), int(betas.shape[1])
    resolved_chunks = chunks or (min(64, max(1, n_samples)), min(4096, max(1, n_loci)))
    array = zarr.create_array(
        path,
        shape=(n_samples, n_loci),
        chunks=resolved_chunks,
        dtype="float32",
        fill_value=np.nan,
    )
    array[:, :] = betas
    return n_samples, n_loci


def open_betas_zarr(path: Path) -> Any:
    return zarr.open_array(path, mode="r")


def read_sample_index(path: Path) -> pd.DataFrame:
    return pd.read_parquet(path)


def read_locus_index(path: Path) -> pd.DataFrame:
    return pd.read_parquet(path)


def write_matrix_manifest(path: Path, manifest: dict[str, Any]) -> None:
    validate_matrix_manifest(manifest)
    write_json(path, manifest)


def source_file_records(paths: list[Path]) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for path in paths:
        resolved = path.resolve()
        if not str(resolved).startswith("/data/"):
            raise ValueError(f"source file must be under /data/: {resolved}")
        records.append(
            {
                "path": str(resolved),
                "sha256": sha256_file(resolved),
                "byte_size": int(resolved.stat().st_size),
            }
        )
    return records
