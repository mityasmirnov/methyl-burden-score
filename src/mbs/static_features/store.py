"""Static feature artifact store (embeddings.zarr + loci.parquet + artifact.json)."""

from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import zarr

from mbs.static_features.manifest import write_static_feature_manifest

LOCI_COLUMNS = (
    "embedding_row",
    "locus_id",
    "canonical_key",
    "source_location_key",
    "source_embedding_row",
    "mapping_status",
)


@dataclass(frozen=True, slots=True)
class StaticFeatureStorePaths:
    root: Path
    embeddings_path: Path
    loci_path: Path
    artifact_path: Path


def static_feature_store_paths(output_dir: Path) -> StaticFeatureStorePaths:
    root = output_dir.resolve()
    return StaticFeatureStorePaths(
        root=root,
        embeddings_path=root / "embeddings.zarr",
        loci_path=root / "loci.parquet",
        artifact_path=root / "artifact.json",
    )


def write_embeddings_zarr(
    path: Path,
    embeddings: np.ndarray,
    *,
    chunks: tuple[int, int] | None = None,
) -> tuple[int, int]:
    """Write ``[n_mapped, dim]`` float16 embeddings to Zarr."""
    if embeddings.ndim != 2:
        raise ValueError(f"embeddings must be 2-D, got shape {embeddings.shape}")
    if embeddings.dtype != np.float16:
        embeddings = embeddings.astype(np.float16, copy=False)
    if path.exists():
        shutil.rmtree(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    n_rows, n_dims = int(embeddings.shape[0]), int(embeddings.shape[1])
    resolved_chunks = chunks or (min(8192, max(1, n_rows)), min(128, max(1, n_dims)))
    array = zarr.create_array(
        path,
        shape=(n_rows, n_dims),
        chunks=resolved_chunks,
        dtype="float16",
    )
    array[:, :] = embeddings
    return n_rows, n_dims


def open_embeddings_zarr(path: Path) -> Any:
    return zarr.open_array(path, mode="r")


def write_loci_index(path: Path, loci: pd.DataFrame) -> pd.DataFrame:
    missing = [col for col in LOCI_COLUMNS if col not in loci.columns]
    if missing:
        raise ValueError(f"loci index missing columns: {missing}")
    frame = loci.loc[:, list(LOCI_COLUMNS)].copy()
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(path, index=False)
    return frame


def read_loci_index(path: Path) -> pd.DataFrame:
    return pd.read_parquet(path)


def write_artifact(path: Path, manifest: dict[str, Any]) -> None:
    write_static_feature_manifest(path, manifest)
