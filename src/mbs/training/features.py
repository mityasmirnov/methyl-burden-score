"""Per-sample methylation + static feature assembly for flat training."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from mbs.training.locus_gene import LocusGeneIndex


@dataclass(frozen=True, slots=True)
class SampleFeatureBundle:
    """Observed CpG→gene edge features for one sample."""

    cpg_features: np.ndarray  # float32 [n_edges_obs, feat_dim]
    cpg_to_gene: np.ndarray  # int64 [n_edges_obs]
    n_observed_edges: int
    n_dropped_nan_beta: int
    n_dropped_no_static: int


def beta_to_m_value(beta: np.ndarray, *, epsilon: float = 0.001) -> np.ndarray:
    """Convert beta values in ``(0, 1)`` to M-values; NaNs preserved."""
    if epsilon <= 0:
        raise ValueError("epsilon must be positive")
    beta = np.asarray(beta, dtype=np.float64)
    clipped = np.clip(beta, epsilon, 1.0 - epsilon)
    m_value = np.log2(clipped / (1.0 - clipped))
    m_value = np.where(np.isfinite(beta), m_value, np.nan)
    return m_value.astype(np.float32, copy=False)


def gather_sample_features(
    *,
    beta_row: np.ndarray,
    static_by_col: np.ndarray,
    static_valid: np.ndarray,
    locus_gene: LocusGeneIndex,
    epsilon: float = 0.001,
    include_m_value: bool = True,
) -> SampleFeatureBundle:
    """Build ragged flat features for one sample.

    ``static_by_col`` is ``[n_study_loci, static_dim]``; ``static_valid`` is
    boolean ``[n_study_loci]``. Edges with non-finite beta or invalid static
    are dropped (no silent zero-fill of missing betas).
    """
    beta_row = np.asarray(beta_row, dtype=np.float32)
    if beta_row.ndim != 1:
        raise ValueError("beta_row must be 1-D")
    if beta_row.shape[0] < locus_gene.n_study_loci:
        raise ValueError(
            f"beta_row length {beta_row.shape[0]} < n_study_loci {locus_gene.n_study_loci}"
        )

    cols = locus_gene.edge_col_index
    genes = locus_gene.edge_gene_index
    betas = beta_row[cols]
    finite = np.isfinite(betas)
    static_ok = static_valid[cols]
    keep = finite & static_ok
    n_dropped_nan = int((~finite).sum())
    n_dropped_static = int((finite & ~static_ok).sum())

    cols_k = cols[keep]
    genes_k = genes[keep]
    betas_k = betas[keep]
    static_k = static_by_col[cols_k]

    parts: list[np.ndarray] = [betas_k.reshape(-1, 1)]
    if include_m_value:
        m_vals = beta_to_m_value(betas_k, epsilon=epsilon).reshape(-1, 1)
        parts.append(m_vals)
    parts.append(static_k.astype(np.float32, copy=False))
    features = np.concatenate(parts, axis=1).astype(np.float32, copy=False)

    return SampleFeatureBundle(
        cpg_features=features,
        cpg_to_gene=genes_k.astype(np.int64, copy=False),
        n_observed_edges=int(features.shape[0]),
        n_dropped_nan_beta=n_dropped_nan,
        n_dropped_no_static=n_dropped_static,
    )


def build_static_column_table(
    *,
    locus_index_locus_ids: np.ndarray,
    static_loci: pd.DataFrame,
    embeddings: Any,
    n_study_loci: int,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Align static embeddings to study columns.

    Returns ``(static_by_col [n_loci, dim], valid [n_loci], dim)``.
    """
    mapped = static_loci.loc[
        static_loci["mapping_status"] == "mapped", ["locus_id", "embedding_row"]
    ].copy()
    emb = np.asarray(embeddings[:], dtype=np.float32)
    if emb.ndim != 2:
        raise ValueError("embeddings must be 2-D")
    dim = int(emb.shape[1])

    study = pd.DataFrame(
        {
            "col_index": np.arange(n_study_loci, dtype=np.int64),
            "locus_id": locus_index_locus_ids[:n_study_loci].astype(np.uint64, copy=False),
        }
    )
    joined = study.merge(mapped, on="locus_id", how="left")
    static_by_col = np.zeros((n_study_loci, dim), dtype=np.float32)
    valid = joined["embedding_row"].notna().to_numpy()
    if valid.any():
        rows = joined.loc[valid, "embedding_row"].astype(int).to_numpy()
        cols = joined.loc[valid, "col_index"].astype(int).to_numpy()
        static_by_col[cols] = emb[rows]
    return static_by_col, valid, dim
