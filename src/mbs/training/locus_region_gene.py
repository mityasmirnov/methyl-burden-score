"""Build CpG→region→gene indices for HierarchicalDeepSet training.

Preserves typed regulatory roles. Gene-unassigned study loci become singleton
``unassigned`` regions under synthetic gene ``__unassigned__``. Illumina
coordinate-unmapped probes are out of scope (never enter matrix columns).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from mbs.annotation.gencode_regions import REGION_TYPES as _GENCODE_REGION_TYPES

UNASSIGNED_GENE_ID = "__unassigned__"
UNASSIGNED_REGION_TYPE = "unassigned"

HIER_REGION_TYPES: tuple[str, ...] = (*_GENCODE_REGION_TYPES, UNASSIGNED_REGION_TYPE)
REGION_TYPE_TO_ID: dict[str, int] = {name: i for i, name in enumerate(HIER_REGION_TYPES)}


@dataclass(frozen=True, slots=True)
class LocusRegionGeneIndex:
    """Study-column aligned locus→region→gene expansion for hierarchy.

    Multi-gene / multi-region loci produce multiple edges (duplicate feature
    rows at train time), matching the flat multi-gene pattern.
    """

    gene_ids: list[str]
    edge_col_index: np.ndarray  # int64 [n_edges]
    edge_region_index: np.ndarray  # int64 [n_edges]
    region_type_id: np.ndarray  # int64 [n_regions]
    region_to_gene: np.ndarray  # int64 [n_regions]
    region_ids: list[str]
    n_study_loci: int
    n_typed_edges: int
    n_unassigned_regions: int

    @property
    def n_genes(self) -> int:
        return len(self.gene_ids)

    @property
    def n_regions(self) -> int:
        return int(self.region_type_id.shape[0])

    @property
    def n_edges(self) -> int:
        return int(self.edge_col_index.shape[0])

    @property
    def unassigned_gene_index(self) -> int | None:
        if UNASSIGNED_GENE_ID not in self.gene_ids:
            return None
        return self.gene_ids.index(UNASSIGNED_GENE_ID)


def build_locus_region_gene_index(
    *,
    locus_index: pd.DataFrame,
    locus_region_edges: pd.DataFrame,
    regions: pd.DataFrame,
    max_loci: int | None = None,
) -> LocusRegionGeneIndex:
    """Join study loci to typed regions; mint singleton orphans for unassigned."""
    required_locus = {"col_index", "locus_id"}
    missing_locus = required_locus - set(locus_index.columns)
    if missing_locus:
        raise ValueError(f"locus_index missing columns: {sorted(missing_locus)}")
    lr_cols = set(locus_region_edges.columns)
    if "locus_id" not in lr_cols or "region_id" not in lr_cols:
        raise ValueError("locus_region_edges requires locus_id and region_id")
    region_cols = set(regions.columns)
    for col in ("region_id", "gene_id", "region_type"):
        if col not in region_cols:
            raise ValueError(f"regions requires {col}")

    study = locus_index.loc[:, ["col_index", "locus_id"]].sort_values("col_index").copy()
    if max_loci is not None:
        if max_loci < 1:
            raise ValueError("max_loci must be >= 1")
        study = study.iloc[:max_loci].copy()
        n_study_loci = int(max_loci)
    else:
        n_study_loci = len(locus_index)

    region_view = regions.loc[:, ["region_id", "gene_id", "region_type"]].copy()
    region_view["region_id"] = region_view["region_id"].astype(str)
    region_view["gene_id"] = region_view["gene_id"].astype(str)
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
            if str(t) not in REGION_TYPE_TO_ID
        }
    )
    if unknown_types:
        raise ValueError(f"unsupported region_type values: {unknown_types}")

    is_typed = (
        merged["region_id"].notna() & merged["gene_id"].notna() & merged["region_type"].notna()
    )
    assigned_loci = set(merged.loc[is_typed, "locus_id"].tolist())
    typed = merged.loc[is_typed, ["col_index", "region_id", "gene_id", "region_type"]].copy()
    orphan_study = study.loc[~study["locus_id"].isin(assigned_loci)].copy()

    frames: list[pd.DataFrame] = []
    if not typed.empty:
        typed = typed.copy()
        typed["region_id"] = typed["region_id"].astype(str)
        typed["gene_id"] = typed["gene_id"].astype(str)
        typed["region_type"] = typed["region_type"].astype(str)
        frames.append(typed.loc[:, ["col_index", "region_id", "gene_id", "region_type"]])

    n_unassigned = 0
    if not orphan_study.empty:
        orphan = pd.DataFrame(
            {
                "col_index": orphan_study["col_index"].to_numpy(dtype=np.int64),
                "region_id": [f"unassigned:{lid}" for lid in orphan_study["locus_id"].tolist()],
                "gene_id": UNASSIGNED_GENE_ID,
                "region_type": UNASSIGNED_REGION_TYPE,
            }
        )
        n_unassigned = len(orphan)
        frames.append(orphan)

    if not frames:
        raise ValueError("no locus→region edges for study loci")

    edges = pd.concat(frames, ignore_index=True)
    # Unique regions in first-seen order.
    region_table = edges.loc[:, ["region_id", "gene_id", "region_type"]].drop_duplicates(
        subset=["region_id"],
        keep="first",
    )
    region_ids = region_table["region_id"].astype(str).tolist()
    region_key_to_idx = {rid: i for i, rid in enumerate(region_ids)}
    region_type_id = np.asarray(
        [REGION_TYPE_TO_ID[str(t)] for t in region_table["region_type"].tolist()],
        dtype=np.int64,
    )
    region_gene_ids = region_table["gene_id"].astype(str).tolist()

    bio_genes = sorted({g for g in region_gene_ids if g != UNASSIGNED_GENE_ID})
    gene_ids = list(bio_genes)
    if UNASSIGNED_GENE_ID in region_gene_ids:
        gene_ids.append(UNASSIGNED_GENE_ID)
    gene_to_idx = {gid: i for i, gid in enumerate(gene_ids)}
    region_to_gene = np.asarray([gene_to_idx[g] for g in region_gene_ids], dtype=np.int64)

    edge_col_index = edges["col_index"].to_numpy(dtype=np.int64, copy=True)
    edge_region_index = np.asarray(
        [region_key_to_idx[str(r)] for r in edges["region_id"].tolist()],
        dtype=np.int64,
    )
    n_typed_edges = int((edges["region_type"] != UNASSIGNED_REGION_TYPE).sum())

    return LocusRegionGeneIndex(
        gene_ids=gene_ids,
        edge_col_index=edge_col_index,
        edge_region_index=edge_region_index,
        region_type_id=region_type_id,
        region_to_gene=region_to_gene,
        region_ids=region_ids,
        n_study_loci=n_study_loci,
        n_typed_edges=n_typed_edges,
        n_unassigned_regions=n_unassigned,
    )
