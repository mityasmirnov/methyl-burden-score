"""Build flat locus→gene indices for FlatDeepSet (skip region encoder)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass(frozen=True, slots=True)
class LocusGeneIndex:
    """Study-column aligned locus→gene edge expansion.

    For each matrix column with ≥1 gene mapping, ``edge_col_index[i]`` is the
    study column and ``edge_gene_index[i]`` is the gene panel index. Multi-gene
    loci produce multiple edges (duplicate feature rows at train time).
    """

    gene_ids: list[str]
    edge_col_index: np.ndarray  # int64 [n_edges]
    edge_gene_index: np.ndarray  # int64 [n_edges]
    n_study_loci: int

    @property
    def n_genes(self) -> int:
        return len(self.gene_ids)

    @property
    def n_edges(self) -> int:
        return int(self.edge_col_index.shape[0])


def build_locus_gene_index(
    *,
    locus_index: pd.DataFrame,
    locus_region_edges: pd.DataFrame,
    regions: pd.DataFrame,
    max_loci: int | None = None,
) -> LocusGeneIndex:
    """Join study loci to genes via ``locus_region_edges`` + ``regions.gene_id``."""
    required_locus = {"col_index", "locus_id"}
    missing_locus = required_locus - set(locus_index.columns)
    if missing_locus:
        raise ValueError(f"locus_index missing columns: {sorted(missing_locus)}")
    lr_cols = set(locus_region_edges.columns)
    if "locus_id" not in lr_cols or "region_id" not in lr_cols:
        raise ValueError("locus_region_edges requires locus_id and region_id")
    if "region_id" not in regions.columns or "gene_id" not in regions.columns:
        raise ValueError("regions requires region_id and gene_id")

    study = locus_index.loc[:, ["col_index", "locus_id"]].sort_values("col_index").copy()
    if max_loci is not None:
        if max_loci < 1:
            raise ValueError("max_loci must be >= 1")
        study = study.iloc[:max_loci].copy()
        n_study_loci = int(max_loci)
    else:
        n_study_loci = len(locus_index)

    merged = study.merge(locus_region_edges[["locus_id", "region_id"]], on="locus_id", how="inner")
    merged = merged.merge(regions[["region_id", "gene_id"]], on="region_id", how="inner")
    if merged.empty:
        raise ValueError("no locus→gene edges for study loci")

    # Distinct (col, gene) edges — a locus in multiple regions of the same gene collapses
    edges = merged.loc[:, ["col_index", "gene_id"]].drop_duplicates()
    gene_ids = sorted(edges["gene_id"].astype(str).unique().tolist())
    gene_to_idx = {gid: i for i, gid in enumerate(gene_ids)}
    edge_col = edges["col_index"].to_numpy(dtype=np.int64, copy=True)
    edge_gene = np.array([gene_to_idx[str(g)] for g in edges["gene_id"].tolist()], dtype=np.int64)
    return LocusGeneIndex(
        gene_ids=gene_ids,
        edge_col_index=edge_col,
        edge_gene_index=edge_gene,
        n_study_loci=n_study_loci,
    )


def load_graph_tables(graph_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load ``locus_region_edges.parquet`` and ``regions.parquet`` from a graph release."""
    root = graph_dir.resolve()
    lr_path = root / "locus_region_edges.parquet"
    regions_path = root / "regions.parquet"
    if not lr_path.is_file():
        raise FileNotFoundError(f"missing {lr_path}")
    if not regions_path.is_file():
        raise FileNotFoundError(f"missing {regions_path}")
    return pd.read_parquet(lr_path), pd.read_parquet(regions_path)
