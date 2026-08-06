"""Validation helpers for static feature exports."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from mbs.static_features.store import LOCI_COLUMNS

_ALLOWED_DTYPES = {"float16", "float32", "float64"}


def validate_embeddings_array(embeddings: np.ndarray, *, output_dimension: int) -> None:
    if embeddings.ndim != 2:
        raise ValueError(f"embeddings must be 2-D, got {embeddings.shape}")
    if embeddings.shape[1] != output_dimension:
        raise ValueError(f"expected output_dimension={output_dimension}, got {embeddings.shape[1]}")
    if str(embeddings.dtype) not in _ALLOWED_DTYPES:
        raise ValueError(f"unexpected embeddings dtype: {embeddings.dtype}")
    as_float = embeddings.astype(np.float64, copy=False)
    if not np.isfinite(as_float).all():
        n_bad = int((~np.isfinite(as_float)).sum())
        raise ValueError(f"embeddings contain {n_bad} non-finite values")


def validate_loci_frame(loci: pd.DataFrame, *, n_mapped: int) -> dict[str, int]:
    missing_cols = [col for col in LOCI_COLUMNS if col not in loci.columns]
    if missing_cols:
        raise ValueError(f"loci frame missing columns: {missing_cols}")
    mapped = loci["mapping_status"].astype(str) == "mapped"
    n_status_mapped = int(mapped.sum())
    if n_status_mapped != n_mapped:
        raise ValueError(
            f"mapping_status mapped count {n_status_mapped} != embeddings rows {n_mapped}"
        )
    mapped_rows = loci.loc[mapped, "embedding_row"]
    if mapped_rows.isna().any():
        raise ValueError("mapped loci must have embedding_row set")
    expected = np.arange(n_mapped, dtype=np.int64)
    got = mapped_rows.astype("int64").to_numpy()
    if not np.array_equal(got, expected):
        raise ValueError("mapped embedding_row values must be contiguous 0..n_mapped-1")
    unmapped = loci.loc[~mapped]
    if unmapped["embedding_row"].notna().any():
        raise ValueError("unmapped loci must have null embedding_row")
    return {
        "n_loci": len(loci),
        "n_mapped": n_status_mapped,
        "n_missing": int((~mapped).sum()),
    }


def embedding_summary_stats(embeddings: np.ndarray) -> dict[str, Any]:
    values = embeddings.astype(np.float64, copy=False)
    norms = np.linalg.norm(values, axis=1)
    dim_var = values.var(axis=0)
    return {
        "n_rows": int(values.shape[0]),
        "n_dims": int(values.shape[1]),
        "norm_mean": float(norms.mean()) if len(norms) else 0.0,
        "norm_std": float(norms.std()) if len(norms) else 0.0,
        "norm_min": float(norms.min()) if len(norms) else 0.0,
        "norm_max": float(norms.max()) if len(norms) else 0.0,
        "dim_var_mean": float(dim_var.mean()) if len(dim_var) else 0.0,
        "dim_var_min": float(dim_var.min()) if len(dim_var) else 0.0,
        "dim_var_max": float(dim_var.max()) if len(dim_var) else 0.0,
        "n_near_zero_norm": int((norms < 1e-6).sum()) if len(norms) else 0,
    }
