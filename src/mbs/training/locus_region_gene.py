"""Build CpG→region→gene indices for HierarchicalDeepSet training.

Preserves typed regulatory roles for annotated loci. Loci / residual probes
without a clean regulatory assignment stay on the residual path — they are
never nearest-gene assigned and never pooled under ``__unassigned__``.

Graph-v2 multi-system mode (``region_systems``) includes RBS/TBS regions as
panel entities (panel id = region_id); gene-only default keeps v0.1 behavior.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from mbs.annotation.gencode_regions import REGION_TYPES as _GENCODE_REGION_TYPES
from mbs.batch import (
    ANNOTATION_STATUS_AMBIGUOUS,
    ANNOTATION_STATUS_MAPPED,
    ANNOTATION_STATUS_MULTI_MAPPED,
    ANNOTATION_STATUS_UNMAPPED,
)
from mbs.matrix.locus_map import is_residual_canonical_key

# Five GENCODE roles — residual loci are not a region type.
HIER_REGION_TYPES: tuple[str, ...] = _GENCODE_REGION_TYPES
RBS_REGION_TYPES: tuple[str, ...] = ("cgi_island", "cgi_north_shore", "cgi_south_shore")
TBS_REGION_TYPES: tuple[str, ...] = ("cpg_tile",)
REGION_TYPE_TO_ID: dict[str, int] = {name: i for i, name in enumerate(HIER_REGION_TYPES)}

RESIDUAL_PANEL_ID = "__residual__"
VALID_REGION_SYSTEMS = frozenset({"gene", "rbs", "tbs"})


def region_type_vocab(region_systems: tuple[str, ...] | list[str]) -> tuple[str, ...]:
    """Region-type embedding vocab for the selected systems."""
    systems = tuple(region_systems)
    types: list[str] = []
    if "gene" in systems:
        types.extend(HIER_REGION_TYPES)
    if "rbs" in systems:
        types.extend(RBS_REGION_TYPES)
    if "tbs" in systems:
        types.extend(TBS_REGION_TYPES)
    if not types:
        raise ValueError(f"region_systems must include gene/rbs/tbs, got {systems}")
    return tuple(types)


def region_systems_from_arm(arm: str | None) -> tuple[str, ...]:
    """Map ``model.arm`` to region systems (default gene-only)."""
    if arm is None or arm == "gene":
        return ("gene",)
    if arm == "rbs":
        return ("rbs",)
    if arm == "tbs":
        return ("tbs",)
    if arm == "direct":
        raise ValueError("direct arm does not use locus→region→gene index")
    raise ValueError(f"unknown arm: {arm}")


@dataclass(frozen=True, slots=True)
class LocusRegionGeneIndex:
    """Study-column aligned locus→region→panel expansion + residual columns.

    Multi-gene / multi-region loci produce multiple typed edges (duplicate
    feature rows at train time). Residual columns never enter typed edges.
    For gene system, panel entities are gene_ids; for rbs/tbs, region_ids.
    """

    gene_ids: list[str]  # panel entity ids (genes and/or RBS/TBS region ids)
    edge_col_index: np.ndarray  # int64 [n_edges] — mapped / multi_mapped only
    edge_region_index: np.ndarray  # int64 [n_edges]
    region_type_id: np.ndarray  # int64 [n_regions]
    region_to_gene: np.ndarray  # int64 [n_regions] — panel index
    region_ids: list[str]
    residual_col_index: np.ndarray  # int64 [n_residual_cols]
    column_annotation_status: np.ndarray  # int8 [n_study_loci]
    n_study_loci: int
    n_typed_edges: int
    region_types: tuple[str, ...] = HIER_REGION_TYPES
    region_systems: tuple[str, ...] = ("gene",)

    @property
    def n_genes(self) -> int:
        return len(self.gene_ids)

    @property
    def panel_ids(self) -> list[str]:
        return self.gene_ids

    @property
    def n_regions(self) -> int:
        return int(self.region_type_id.shape[0])

    @property
    def n_edges(self) -> int:
        return int(self.edge_col_index.shape[0])

    @property
    def n_residual_cols(self) -> int:
        return int(self.residual_col_index.shape[0])

    @property
    def n_panel(self) -> int:
        """Panel size plus one residual score slot."""
        return self.n_genes + 1

    @property
    def residual_panel_index(self) -> int:
        return self.n_genes


def build_locus_region_gene_index(
    *,
    locus_index: pd.DataFrame,
    locus_region_edges: pd.DataFrame,
    regions: pd.DataFrame,
    max_loci: int | None = None,
    region_systems: tuple[str, ...] | list[str] = ("gene",),
) -> LocusRegionGeneIndex:
    """Join study loci to typed regions; route orphans to the residual path."""
    systems = tuple(region_systems)
    unknown_sys = set(systems) - VALID_REGION_SYSTEMS
    if unknown_sys:
        raise ValueError(f"unsupported region_systems: {sorted(unknown_sys)}")
    if not systems:
        raise ValueError("region_systems must be non-empty")
    type_vocab = region_type_vocab(systems)
    type_to_id = {name: i for i, name in enumerate(type_vocab)}

    required_locus = {"col_index", "locus_id"}
    missing_locus = required_locus - set(locus_index.columns)
    if missing_locus:
        raise ValueError(f"locus_index missing columns: {sorted(missing_locus)}")
    lr_cols = set(locus_region_edges.columns)
    if "locus_id" not in lr_cols or "region_id" not in lr_cols:
        raise ValueError("locus_region_edges requires locus_id and region_id")
    region_cols = set(regions.columns)
    for col in ("region_id", "region_type"):
        if col not in region_cols:
            raise ValueError(f"regions requires {col}")
    if "gene_id" not in region_cols and "gene" in systems:
        raise ValueError("regions requires gene_id for gene system")

    selected = regions.copy()
    if "region_system" in region_cols:
        selected["region_system"] = selected["region_system"].fillna("gene")
    else:
        selected["region_system"] = "gene"
    selected = selected.loc[selected["region_system"].isin(systems)].copy()
    if selected.empty:
        raise ValueError(f"no regions for region_systems={list(systems)}")

    selected = selected.reset_index(drop=True)
    # Panel entity: gene_id for gene; region_id for rbs/tbs.
    panel_ids: list[str] = []
    keep: list[bool] = []
    for rec in selected.itertuples(index=False):
        sys = str(rec.region_system)
        if sys == "gene":
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
        raise ValueError(f"no typed regions after panel filter for systems={list(systems)}")

    study = locus_index.loc[:, ["col_index", "locus_id"]].sort_values("col_index").copy()
    if "canonical_key" in locus_index.columns:
        study = study.merge(
            locus_index.loc[:, ["col_index", "canonical_key"]],
            on="col_index",
            how="left",
        )
    else:
        study["canonical_key"] = ""
    if max_loci is not None:
        if max_loci < 1:
            raise ValueError("max_loci must be >= 1")
        study = study.iloc[:max_loci].copy()
        n_study_loci = int(max_loci)
    else:
        n_study_loci = len(locus_index)

    region_view = selected.loc[:, ["region_id", "panel_id", "region_type"]].copy()
    region_view["region_id"] = region_view["region_id"].astype(str)
    region_view["panel_id"] = region_view["panel_id"].astype(str)
    region_view["region_type"] = region_view["region_type"].astype(str)

    lr = locus_region_edges.loc[:, ["locus_id", "region_id"]].copy()
    lr["region_id"] = lr["region_id"].astype(str)

    merged = study.merge(lr, on="locus_id", how="left")
    merged = merged.merge(region_view, on="region_id", how="left")
    if merged.empty:
        raise ValueError("no study loci after join")

    unknown_types = sorted(
        {
            str(t)
            for t in merged["region_type"].dropna().unique().tolist()
            if str(t) not in type_to_id
        }
    )
    if unknown_types:
        raise ValueError(f"unsupported region_type values: {unknown_types}")

    is_typed = (
        merged["region_id"].notna() & merged["panel_id"].notna() & merged["region_type"].notna()
    )
    is_ambiguous_edge = merged["region_id"].notna() & ~is_typed

    typed = merged.loc[is_typed, ["col_index", "region_id", "panel_id", "region_type"]].copy()
    typed_cols = set(typed["col_index"].astype(int).tolist()) if not typed.empty else set()
    ambiguous_cols = (
        set(merged.loc[is_ambiguous_edge, "col_index"].astype(int).tolist()) - typed_cols
    )

    col_panel_n = np.zeros(n_study_loci, dtype=np.int32)
    if not typed.empty:
        panel_counts = typed.groupby("col_index")["panel_id"].nunique().astype(int)
        idx = panel_counts.index.to_numpy(dtype=np.int64)
        vals = panel_counts.to_numpy(dtype=np.int32)
        keep = (idx >= 0) & (idx < n_study_loci)
        col_panel_n[idx[keep]] = vals[keep]

    status = np.full(n_study_loci, ANNOTATION_STATUS_UNMAPPED, dtype=np.int8)
    cols = study["col_index"].to_numpy(dtype=np.int64, copy=False)
    keys = study["canonical_key"].astype(str).fillna("").to_numpy()
    in_range = (cols >= 0) & (cols < n_study_loci)
    cols_ok = cols[in_range]
    keys_ok = keys[in_range]

    residual_key = np.fromiter(
        (is_residual_canonical_key(str(k)) for k in keys_ok),
        dtype=np.bool_,
        count=int(keys_ok.shape[0]),
    )
    residual_col_mask = np.zeros(n_study_loci, dtype=np.bool_)
    residual_col_mask[cols_ok[residual_key]] = True

    ambiguous_mask = np.zeros(n_study_loci, dtype=np.bool_)
    if ambiguous_cols:
        amb = np.fromiter(
            (c for c in ambiguous_cols if 0 <= c < n_study_loci and c not in typed_cols),
            dtype=np.int64,
        )
        if amb.size:
            ambiguous_mask[amb] = True

    mapped_mask = (col_panel_n == 1) & ~residual_col_mask & ~ambiguous_mask
    multi_mask = (col_panel_n >= 2) & ~residual_col_mask & ~ambiguous_mask
    status[mapped_mask] = ANNOTATION_STATUS_MAPPED
    status[multi_mask] = ANNOTATION_STATUS_MULTI_MAPPED
    status[ambiguous_mask] = ANNOTATION_STATUS_AMBIGUOUS
    status[residual_col_mask & ~ambiguous_mask] = ANNOTATION_STATUS_UNMAPPED

    residual_path_mask = residual_col_mask | ambiguous_mask | (col_panel_n == 0)
    residual_col_index = np.flatnonzero(residual_path_mask).astype(np.int64)

    if typed.empty and residual_col_index.size == 0:
        raise ValueError("no locus→region edges and no residual columns for study loci")

    if typed.empty:
        return LocusRegionGeneIndex(
            gene_ids=[],
            edge_col_index=np.zeros(0, dtype=np.int64),
            edge_region_index=np.zeros(0, dtype=np.int64),
            region_type_id=np.zeros(0, dtype=np.int64),
            region_to_gene=np.zeros(0, dtype=np.int64),
            region_ids=[],
            residual_col_index=residual_col_index,
            column_annotation_status=status,
            n_study_loci=n_study_loci,
            n_typed_edges=0,
            region_types=type_vocab,
            region_systems=systems,
        )

    typed["region_id"] = typed["region_id"].astype(str)
    typed["panel_id"] = typed["panel_id"].astype(str)
    typed["region_type"] = typed["region_type"].astype(str)
    edges = typed.loc[:, ["col_index", "region_id", "panel_id", "region_type"]]

    region_table = edges.loc[:, ["region_id", "panel_id", "region_type"]].drop_duplicates(
        subset=["region_id"],
        keep="first",
    )
    region_ids = region_table["region_id"].astype(str).tolist()
    region_key_to_idx = {rid: i for i, rid in enumerate(region_ids)}
    region_type_id = np.asarray(
        [type_to_id[str(t)] for t in region_table["region_type"].tolist()],
        dtype=np.int64,
    )
    region_panel_ids = region_table["panel_id"].astype(str).tolist()
    gene_ids = sorted(set(region_panel_ids))
    gene_to_idx = {gid: i for i, gid in enumerate(gene_ids)}
    region_to_gene = np.asarray([gene_to_idx[g] for g in region_panel_ids], dtype=np.int64)

    edge_col_index = edges["col_index"].to_numpy(dtype=np.int64, copy=True)
    edge_region_index = np.asarray(
        [region_key_to_idx[str(r)] for r in edges["region_id"].tolist()],
        dtype=np.int64,
    )

    return LocusRegionGeneIndex(
        gene_ids=gene_ids,
        edge_col_index=edge_col_index,
        edge_region_index=edge_region_index,
        region_type_id=region_type_id,
        region_to_gene=region_to_gene,
        region_ids=region_ids,
        residual_col_index=residual_col_index,
        column_annotation_status=status,
        n_study_loci=n_study_loci,
        n_typed_edges=int(edge_col_index.shape[0]),
        region_types=type_vocab,
        region_systems=systems,
    )


def locus_region_gene_col_filter(
    locus_region: LocusRegionGeneIndex,
    cols: np.ndarray,
) -> LocusRegionGeneIndex:
    """Keep typed edges on ``cols`` only; drop residual path (7G′ gene_cols parity)."""
    allowed = frozenset(int(c) for c in np.asarray(cols, dtype=np.int64).tolist())
    if not allowed:
        raise ValueError("cols must be non-empty")
    mask = np.fromiter(
        (int(c) in allowed for c in locus_region.edge_col_index),
        dtype=bool,
        count=locus_region.n_edges,
    )
    edge_col = locus_region.edge_col_index[mask]
    edge_reg = locus_region.edge_region_index[mask]
    return LocusRegionGeneIndex(
        gene_ids=locus_region.gene_ids,
        edge_col_index=edge_col,
        edge_region_index=edge_reg,
        region_type_id=locus_region.region_type_id,
        region_to_gene=locus_region.region_to_gene,
        region_ids=locus_region.region_ids,
        residual_col_index=np.zeros(0, dtype=np.int64),
        column_annotation_status=locus_region.column_annotation_status,
        n_study_loci=locus_region.n_study_loci,
        n_typed_edges=int(edge_col.shape[0]),
        region_types=locus_region.region_types,
        region_systems=locus_region.region_systems,
    )
