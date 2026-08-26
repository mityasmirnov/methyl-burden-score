"""Milestone 7F cascade assignment: typed RBS → gene / orphan; leftover → direct.

Ignores ``region_system=tbs``. Nearest-gene allocates typed RBS with null
``gene_id`` onto a gene (MBS); leftover CpGs are never nearest-gene collapsed
(ADR 0004 / ADR 0009).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from mbs.training.locus_region_gene import (
    RBS_REGION_TYPES,
    HIER_REGION_TYPES,
    region_type_vocab,
)

TYPED_SYSTEMS: tuple[str, ...] = ("gene", "rbs")
ORPHAN_GENE_INDEX = -1


@dataclass(frozen=True, slots=True)
class CascadeAssignment:
    """CpG→typed-region index with gene allocation and leftover-direct columns."""

    gene_ids: list[str]
    region_ids: list[str]
    region_type_id: np.ndarray  # int64 [n_regions]
    region_to_gene: np.ndarray  # int64 [n_regions]; ORPHAN_GENE_INDEX if orphan
    orphan_region_mask: np.ndarray  # bool [n_regions]
    edge_col_index: np.ndarray  # int64 [n_edges]
    edge_region_index: np.ndarray  # int64 [n_edges]
    direct_col_index: np.ndarray  # int64 [n_direct]
    region_types: tuple[str, ...]
    n_study_loci: int
    allocated_gene_id: list[str | None]  # parallel to region_ids

    @property
    def n_genes(self) -> int:
        return len(self.gene_ids)

    @property
    def n_regions(self) -> int:
        return len(self.region_ids)

    @property
    def n_orphan_rbs(self) -> int:
        return int(self.orphan_region_mask.sum())

    @property
    def n_direct(self) -> int:
        return int(self.direct_col_index.shape[0])

    @property
    def orphan_region_indices(self) -> np.ndarray:
        return np.flatnonzero(self.orphan_region_mask).astype(np.int64)

    @property
    def orphan_region_ids(self) -> list[str]:
        return [self.region_ids[i] for i in self.orphan_region_indices.tolist()]


def nearest_gene_on_chromosome(
    chromosome: str,
    midpoint: int,
    genes: pd.DataFrame,
) -> str | None:
    """Return nearest gene_id on ``chromosome`` by distance to gene interval; else None."""
    if genes.empty or "chromosome" not in genes.columns:
        return None
    chrom = str(chromosome)
    g = genes.loc[genes["chromosome"].astype(str) == chrom]
    if g.empty:
        return None
    starts = g["start"].to_numpy(dtype=np.int64, copy=False)
    ends = g["end"].to_numpy(dtype=np.int64, copy=False)
    # Distance 0 if midpoint inside [start, end]; else to nearest endpoint.
    dist = np.where(
        (midpoint >= starts) & (midpoint <= ends),
        0,
        np.minimum(np.abs(midpoint - starts), np.abs(midpoint - ends)),
    )
    best = int(np.argmin(dist))
    return str(g.iloc[best]["gene_id"])


def allocate_rbs_genes(
    regions: pd.DataFrame,
    genes: pd.DataFrame,
) -> pd.Series:
    """Map each region row to allocated gene_id (typed gene_id or nearest-gene).

    Regions already carrying a non-null ``gene_id`` keep it. Null-gene typed
    RBS regions receive same-chromosome nearest gene; failure → None (orphan).
    """
    out: list[str | None] = []
    for rec in regions.itertuples(index=False):
        gid = getattr(rec, "gene_id", None)
        if gid is not None and not (isinstance(gid, float) and np.isnan(gid)) and str(gid) not in (
            "",
            ".",
            "None",
            "nan",
        ):
            out.append(str(gid))
            continue
        chrom = str(getattr(rec, "chromosome", ""))
        start = int(getattr(rec, "start", 0) or 0)
        end = int(getattr(rec, "end", start) or start)
        mid = (start + end) // 2
        out.append(nearest_gene_on_chromosome(chrom, mid, genes))
    return pd.Series(out, index=regions.index, dtype=object)


def build_cascade_assignment(
    *,
    locus_index: pd.DataFrame,
    locus_region_edges: pd.DataFrame,
    regions: pd.DataFrame,
    genes: pd.DataFrame,
    max_loci: int | None = None,
) -> CascadeAssignment:
    """Build 7F assignment: gene+rbs only; leftover (incl. former TBS) → direct."""
    required_locus = {"col_index", "locus_id"}
    missing_locus = required_locus - set(locus_index.columns)
    if missing_locus:
        raise ValueError(f"locus_index missing columns: {sorted(missing_locus)}")
    if "locus_id" not in locus_region_edges.columns or "region_id" not in locus_region_edges.columns:
        raise ValueError("locus_region_edges requires locus_id and region_id")
    for col in ("region_id", "region_type"):
        if col not in regions.columns:
            raise ValueError(f"regions requires {col}")

    selected = regions.copy()
    if "region_system" in selected.columns:
        selected["region_system"] = selected["region_system"].fillna("gene")
    else:
        # Infer: GENCODE roles → gene; CGI types → rbs; tiles → tbs (dropped).
        gene_types = set(HIER_REGION_TYPES)
        rbs_types = set(RBS_REGION_TYPES)

        def _infer_system(rt: object) -> str:
            s = str(rt)
            if s in gene_types:
                return "gene"
            if s in rbs_types:
                return "rbs"
            if s == "cpg_tile":
                return "tbs"
            return "rbs"

        selected["region_system"] = selected["region_type"].map(_infer_system)

    # Drop TBS entirely for scoring.
    selected = selected.loc[selected["region_system"].isin(TYPED_SYSTEMS)].copy()
    if selected.empty:
        raise ValueError("no gene/rbs regions for cascade assignment")

    selected["allocated_gene_id"] = allocate_rbs_genes(selected, genes)

    type_vocab = region_type_vocab(TYPED_SYSTEMS)
    type_to_id = {name: i for i, name in enumerate(type_vocab)}

    study = locus_index.loc[:, ["col_index", "locus_id"]].sort_values("col_index").copy()
    if max_loci is not None:
        if max_loci < 1:
            raise ValueError("max_loci must be >= 1")
        study = study.iloc[:max_loci].copy()
        n_study_loci = int(max_loci)
    else:
        n_study_loci = int(len(locus_index))

    region_view = selected.loc[
        :,
        [c for c in ("region_id", "region_type", "allocated_gene_id", "region_system") if c in selected.columns],
    ].copy()
    region_view["region_id"] = region_view["region_id"].astype(str)
    region_view["region_type"] = region_view["region_type"].astype(str)

    lr = locus_region_edges.loc[:, ["locus_id", "region_id"]].copy()
    lr["region_id"] = lr["region_id"].astype(str)
    # Drop edges that point at TBS / unknown regions (left join miss).
    merged = study.merge(lr, on="locus_id", how="left")
    merged = merged.merge(region_view, on="region_id", how="left")

    is_typed = merged["region_id"].notna() & merged["region_type"].notna()
    typed = merged.loc[is_typed, ["col_index", "region_id", "region_type", "allocated_gene_id"]].copy()
    typed_cols = set(typed["col_index"].astype(int).tolist()) if not typed.empty else set()

    # Direct = study columns with no typed gene/rbs edge (includes former TBS-only).
    all_cols = study["col_index"].to_numpy(dtype=np.int64, copy=False)
    in_range = (all_cols >= 0) & (all_cols < n_study_loci)
    study_cols = set(int(c) for c in all_cols[in_range].tolist())
    direct_cols = sorted(study_cols - typed_cols)
    direct_col_index = np.asarray(direct_cols, dtype=np.int64)

    if typed.empty:
        return CascadeAssignment(
            gene_ids=[],
            region_ids=[],
            region_type_id=np.zeros(0, dtype=np.int64),
            region_to_gene=np.zeros(0, dtype=np.int64),
            orphan_region_mask=np.zeros(0, dtype=bool),
            edge_col_index=np.zeros(0, dtype=np.int64),
            edge_region_index=np.zeros(0, dtype=np.int64),
            direct_col_index=direct_col_index,
            region_types=type_vocab,
            n_study_loci=n_study_loci,
            allocated_gene_id=[],
        )

    unknown_types = sorted(
        {str(t) for t in typed["region_type"].unique().tolist() if str(t) not in type_to_id}
    )
    if unknown_types:
        raise ValueError(f"unsupported region_type values: {unknown_types}")

    region_table = (
        typed.loc[:, ["region_id", "region_type", "allocated_gene_id"]]
        .drop_duplicates(subset=["region_id"], keep="first")
        .reset_index(drop=True)
    )
    region_ids = region_table["region_id"].astype(str).tolist()
    region_key_to_idx = {rid: i for i, rid in enumerate(region_ids)}
    region_type_id = np.asarray(
        [type_to_id[str(t)] for t in region_table["region_type"].tolist()],
        dtype=np.int64,
    )
    allocated = [
        None if (g is None or (isinstance(g, float) and np.isnan(g)) or str(g) in ("", "None", "nan"))
        else str(g)
        for g in region_table["allocated_gene_id"].tolist()
    ]
    gene_ids = sorted({g for g in allocated if g is not None})
    gene_to_idx = {gid: i for i, gid in enumerate(gene_ids)}
    region_to_gene = np.asarray(
        [gene_to_idx[g] if g is not None else ORPHAN_GENE_INDEX for g in allocated],
        dtype=np.int64,
    )
    orphan_region_mask = region_to_gene < 0

    edge_col_index = typed["col_index"].to_numpy(dtype=np.int64, copy=True)
    edge_region_index = np.asarray(
        [region_key_to_idx[str(r)] for r in typed["region_id"].tolist()],
        dtype=np.int64,
    )

    return CascadeAssignment(
        gene_ids=gene_ids,
        region_ids=region_ids,
        region_type_id=region_type_id,
        region_to_gene=region_to_gene,
        orphan_region_mask=orphan_region_mask,
        edge_col_index=edge_col_index,
        edge_region_index=edge_region_index,
        direct_col_index=direct_col_index,
        region_types=type_vocab,
        n_study_loci=n_study_loci,
        allocated_gene_id=allocated,
    )
