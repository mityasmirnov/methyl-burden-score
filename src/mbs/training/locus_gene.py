"""Build flat locus→gene (or RBS/TBS panel) indices for FlatDeepSet."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from mbs.training.locus_region_gene import VALID_REGION_SYSTEMS, region_systems_from_arm

__all__ = (
    "LocusGeneIndex",
    "build_locus_gene_index",
    "locus_gene_col_filter",
    "load_graph_tables",
    "region_systems_from_arm",
)


@dataclass(frozen=True, slots=True)
class LocusGeneIndex:
    """Study-column aligned locus→panel edge expansion.

    For each matrix column with ≥1 mapping, ``edge_col_index[i]`` is the
    study column and ``edge_gene_index[i]`` is the panel index. Multi-entity
    loci produce multiple edges (duplicate feature rows at train time).
    Panel entities are genes (gene system) or region ids (rbs/tbs).
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
    region_systems: tuple[str, ...] | list[str] = ("gene",),
) -> LocusGeneIndex:
    """Join study loci to panel entities via ``locus_region_edges`` + ``regions``.

    Gene system uses ``gene_id``; rbs/tbs use ``region_id`` as the panel id.
    """
    systems = tuple(region_systems)
    unknown = set(systems) - VALID_REGION_SYSTEMS
    if unknown:
        raise ValueError(f"unsupported region_systems: {sorted(unknown)}")
    if not systems:
        raise ValueError("region_systems must be non-empty")

    required_locus = {"col_index", "locus_id"}
    missing_locus = required_locus - set(locus_index.columns)
    if missing_locus:
        raise ValueError(f"locus_index missing columns: {sorted(missing_locus)}")
    lr_cols = set(locus_region_edges.columns)
    if "locus_id" not in lr_cols or "region_id" not in lr_cols:
        raise ValueError("locus_region_edges requires locus_id and region_id")
    if "region_id" not in regions.columns:
        raise ValueError("regions requires region_id")

    study = locus_index.loc[:, ["col_index", "locus_id"]].sort_values("col_index").copy()
    if max_loci is not None:
        if max_loci < 1:
            raise ValueError("max_loci must be >= 1")
        study = study.iloc[:max_loci].copy()
        n_study_loci = int(max_loci)
    else:
        n_study_loci = len(locus_index)

    selected = regions.copy()
    if "region_system" in selected.columns:
        selected["region_system"] = selected["region_system"].fillna("gene")
    else:
        selected["region_system"] = "gene"
    selected = selected.loc[selected["region_system"].isin(systems)].copy()
    if selected.empty:
        raise ValueError(f"no regions for region_systems={list(systems)}")

    selected = selected.reset_index(drop=True)
    panel_ids: list[str] = []
    keep: list[bool] = []
    for rec in selected.itertuples(index=False):
        sys = str(rec.region_system)
        if sys == "gene":
            if "gene_id" not in selected.columns:
                raise ValueError("regions requires gene_id for gene system")
            gid = rec.gene_id
            if gid is None or pd.isna(gid):
                keep.append(False)
                panel_ids.append("")
            else:
                keep.append(True)
                panel_ids.append(str(gid))
        else:
            keep.append(True)
            panel_ids.append(str(rec.region_id))
    keep_arr = np.asarray(keep, dtype=bool)
    selected = selected.loc[keep_arr].copy()
    selected["panel_id"] = [p for p, k in zip(panel_ids, keep, strict=True) if k]
    if selected.empty:
        raise ValueError(f"no panel entities for region_systems={list(systems)}")

    merged = study.merge(locus_region_edges[["locus_id", "region_id"]], on="locus_id", how="inner")
    merged = merged.merge(
        selected[["region_id", "panel_id"]],
        on="region_id",
        how="inner",
    )
    if merged.empty:
        raise ValueError(f"no locus→panel edges for region_systems={list(systems)}")

    edges = merged.loc[:, ["col_index", "panel_id"]].drop_duplicates()
    gene_ids = sorted(edges["panel_id"].astype(str).unique().tolist())
    gene_to_idx = {gid: i for i, gid in enumerate(gene_ids)}
    edge_col = edges["col_index"].to_numpy(dtype=np.int64, copy=True)
    edge_gene = np.array([gene_to_idx[str(g)] for g in edges["panel_id"].tolist()], dtype=np.int64)
    return LocusGeneIndex(
        gene_ids=gene_ids,
        edge_col_index=edge_col,
        edge_gene_index=edge_gene,
        n_study_loci=n_study_loci,
    )


def locus_gene_col_filter(locus_gene: LocusGeneIndex, cols: np.ndarray) -> LocusGeneIndex:
    """Keep only edges whose study column is in ``cols`` (7G′ gene_cols parity)."""
    allowed = frozenset(int(c) for c in np.asarray(cols, dtype=np.int64).tolist())
    if not allowed:
        raise ValueError("cols must be non-empty")
    mask = np.fromiter(
        (int(c) in allowed for c in locus_gene.edge_col_index),
        dtype=bool,
        count=locus_gene.n_edges,
    )
    return LocusGeneIndex(
        gene_ids=locus_gene.gene_ids,
        edge_col_index=locus_gene.edge_col_index[mask],
        edge_gene_index=locus_gene.edge_gene_index[mask],
        n_study_loci=locus_gene.n_study_loci,
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
