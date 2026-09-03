"""Write 7F score matrices (MBS / orphan RBS / direct); no TBS."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pandas as pd
import zarr

from mbs.scoring.orientation import score_manifest

FusionBlockMode = Literal["full", "mbs_direct", "mbs_only"]


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
    direct_cpg: np.ndarray | None = None,
    direct_locus_ids: list[str] | None = None,
    all_gene_rbs: np.ndarray | None = None,
    all_gene_rbs_present: np.ndarray | None = None,
    all_gene_region_ids: list[str] | None = None,
    all_gene_region_gene_ids: list[str | None] | None = None,
    all_gene_region_types: list[str] | None = None,
    allocation_policy: str | None = None,
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

    if direct_cpg is not None:
        dcpg_a = np.asarray(direct_cpg, dtype=np.float32)
        if direct_locus_ids is None:
            raise ValueError("direct_locus_ids required when direct_cpg is set")
        if dcpg_a.shape != (n, len(direct_locus_ids)):
            raise ValueError(
                f"direct_cpg shape {dcpg_a.shape} != ({n}, {len(direct_locus_ids)})"
            )
        _write_array(out_dir / "direct_cpg.zarr", dcpg_a)
        pd.DataFrame(
            {
                "locus_id": direct_locus_ids,
                "col_index": np.arange(len(direct_locus_ids), dtype=np.int64),
            }
        ).to_parquet(out_dir / "direct_locus_index.parquet", index=False)

    if all_gene_rbs is not None:
        agr = np.asarray(all_gene_rbs, dtype=np.float32)
        agr_present = (
            np.asarray(all_gene_rbs_present, dtype=bool)
            if all_gene_rbs_present is not None
            else np.ones(agr.shape, dtype=bool)
        )
        region_ids = list(all_gene_region_ids or [])
        if agr.shape != (n, len(region_ids)):
            raise ValueError(
                f"all_gene_rbs shape {agr.shape} != ({n}, {len(region_ids)})"
            )
        if agr_present.shape != agr.shape:
            raise ValueError("all_gene_rbs_present shape mismatch")
        gene_ids_col = list(all_gene_region_gene_ids or [None] * len(region_ids))
        type_col = list(all_gene_region_types or ["unknown"] * len(region_ids))
        if len(gene_ids_col) != len(region_ids) or len(type_col) != len(region_ids):
            raise ValueError("all_gene region index length mismatch")
        _write_array(out_dir / "all_gene_rbs.zarr", agr)
        _write_array(out_dir / "all_gene_rbs_present.zarr", agr_present.astype(np.uint8))
        pd.DataFrame(
            {
                "region_id": region_ids,
                "gene_id": gene_ids_col,
                "region_type": type_col,
                "column_index": np.arange(len(region_ids), dtype=np.int64),
                "allocation_policy": allocation_policy or "unknown",
            }
        ).to_parquet(out_dir / "all_gene_region_index.parquet", index=False)

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
    if direct_cpg is not None:
        manifest["n_direct_loci"] = len(direct_locus_ids or [])
        manifest["direct_cpg"] = True
    if all_gene_rbs is not None:
        manifest["n_all_gene_rbs"] = int(np.asarray(all_gene_rbs).shape[1])
        manifest["all_gene_rbs"] = True
        if allocation_policy is not None:
            manifest["allocation_policy"] = allocation_policy
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
    out: dict[str, np.ndarray] = {"mbs": mbs, "orphan_rbs": rbs, "direct": direct}
    all_gene_path = root / "all_gene_rbs.zarr"
    if all_gene_path.is_dir() or all_gene_path.is_file():
        out["all_gene_rbs"] = np.asarray(
            zarr.open_array(str(all_gene_path), mode="r"), dtype=np.float32
        )
        present_path = root / "all_gene_rbs_present.zarr"
        if present_path.is_dir() or present_path.is_file():
            out["all_gene_rbs_present"] = np.asarray(
                zarr.open_array(str(present_path), mode="r"), dtype=bool
            )
    return out


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
