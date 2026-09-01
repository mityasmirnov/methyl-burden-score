"""Write 7F score matrices (MBS / orphan RBS / direct); no TBS."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Literal

FusionBlockMode = Literal["full", "mbs_direct", "mbs_only"]

import numpy as np
import pandas as pd
import zarr

from mbs.scoring.orientation import score_manifest


def _write_array(path: Path, data: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        import shutil

        shutil.rmtree(path)
    arr = zarr.create_array(
        path,
        shape=data.shape,
        chunks=data.shape,
        dtype=data.dtype,
        compressors="auto",
        overwrite=True,
    )
    arr[...] = data


def write_cascade_score_dir(
    out_dir: Path,
    *,
    sample_ids: list[str],
    gene_ids: list[str],
    orphan_region_ids: list[str],
    mbs: np.ndarray,
    gene_present: np.ndarray,
    orphan_rbs: np.ndarray,
    direct_contrib: np.ndarray,
    direct_task_names: list[str],
    score_polarity: str = "hyper_aligned",
    fold_id: str | None = None,
    restart_id: str | None = None,
    extra_manifest: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Persist sample×score Zarrs under ``out_dir`` (DATA_CONTRACT 7F layout)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    mbs_a = np.asarray(mbs, dtype=np.float32)
    present_a = np.asarray(gene_present, dtype=bool)
    rbs_a = np.asarray(orphan_rbs, dtype=np.float32)
    direct_a = np.asarray(direct_contrib, dtype=np.float32)
    n = len(sample_ids)
    if mbs_a.shape != (n, len(gene_ids)):
        raise ValueError(f"mbs shape {mbs_a.shape} != ({n}, {len(gene_ids)})")
    if present_a.shape != mbs_a.shape:
        raise ValueError("gene_present shape mismatch")
    if rbs_a.shape != (n, len(orphan_region_ids)):
        raise ValueError(f"orphan_rbs shape {rbs_a.shape} != ({n}, {len(orphan_region_ids)})")
    if direct_a.ndim != 2 or direct_a.shape[0] != n:
        raise ValueError("direct_contrib must be [n_samples, n_tasks]")
    if direct_a.shape[1] != len(direct_task_names):
        raise ValueError("direct_task_names length mismatch")

    _write_array(out_dir / "mbs.zarr", mbs_a)
    _write_array(out_dir / "gene_present.zarr", present_a.astype(np.uint8))
    _write_array(out_dir / "rbs.zarr", rbs_a)
    _write_array(out_dir / "direct_contrib.zarr", direct_a)

    pd.DataFrame({"sample_id": sample_ids, "row_index": np.arange(n, dtype=np.int64)}).to_parquet(
        out_dir / "sample_index.parquet", index=False
    )
    pd.DataFrame({"gene_id": gene_ids, "col_index": np.arange(len(gene_ids), dtype=np.int64)}).to_parquet(
        out_dir / "gene_index.parquet", index=False
    )
    pd.DataFrame(
        {
            "region_id": orphan_region_ids,
            "col_index": np.arange(len(orphan_region_ids), dtype=np.int64),
        }
    ).to_parquet(out_dir / "region_index.parquet", index=False)
    pd.DataFrame({"task": direct_task_names, "col_index": np.arange(len(direct_task_names))}).to_parquet(
        out_dir / "direct_task_index.parquet", index=False
    )

    manifest = score_manifest(
        score_polarity=score_polarity,
        fold_id=fold_id,
        restart_id=restart_id,
    )
    manifest["topology"] = "rbs_gene_direct_7f"
    manifest["n_genes"] = len(gene_ids)
    manifest["n_orphan_rbs"] = len(orphan_region_ids)
    manifest["n_direct_tasks"] = len(direct_task_names)
    manifest["tbs_arm"] = False
    if extra_manifest:
        manifest.update(extra_manifest)
    (out_dir / "score_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest


def load_cascade_score_blocks(score_dir: Path) -> dict[str, np.ndarray]:
    """Load orphan RBS, MBS, and direct blocks for late fusion (no TBS)."""
    root = Path(score_dir)
    mbs = np.asarray(zarr.open_array(str(root / "mbs.zarr"), mode="r"), dtype=np.float32)
    rbs = np.asarray(zarr.open_array(str(root / "rbs.zarr"), mode="r"), dtype=np.float32)
    direct = np.asarray(
        zarr.open_array(str(root / "direct_contrib.zarr"), mode="r"), dtype=np.float32
    )
    return {"mbs": mbs, "orphan_rbs": rbs, "direct": direct}


def fusion_feature_matrix(
    blocks: dict[str, np.ndarray],
    *,
    mode: FusionBlockMode = "full",
) -> np.ndarray:
    """Column-stack score blocks for late fusion (7F layout; no TBS)."""
    if "tbs" in blocks:
        raise ValueError("TBS arm is forbidden in 7F fusion matrix")
    if mode == "full":
        parts = [blocks["orphan_rbs"], blocks["mbs"], blocks["direct"]]
    elif mode == "mbs_direct":
        parts = [blocks["mbs"], blocks["direct"]]
    elif mode == "mbs_only":
        parts = [blocks["mbs"]]
    else:
        raise ValueError(f"unsupported fusion block mode: {mode!r}")
    return np.concatenate([np.asarray(p, dtype=np.float32) for p in parts], axis=1)
