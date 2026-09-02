"""Milestone 7F cascade assignment: typed RBS → gene / orphan; leftover → direct.

Ignores ``region_system=tbs``. Nearest-gene allocates typed RBS with null
``gene_id`` onto a gene (MBS); leftover CpGs are never nearest-gene collapsed
(ADR 0004 / ADR 0009).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
import pandas as pd

from mbs.training.locus_region_gene import (
    RBS_REGION_TYPES,
    HIER_REGION_TYPES,
    region_type_vocab,
)

TYPED_SYSTEMS: tuple[str, ...] = ("gene", "rbs")
ORPHAN_GENE_INDEX = -1
GeneAllocationPolicy = Literal["explicit_only", "bounded_nearest", "legacy_nearest"]
DEFAULT_BOUNDED_NEAREST_BP = 1_000_000


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


def gene_linked_col_index(assignment: CascadeAssignment) -> np.ndarray:
    """Unique CpG column indices on edges allocated to a gene (MBS path)."""
    if assignment.edge_col_index.size == 0:
        return np.zeros(0, dtype=np.int64)
    gene_edge = assignment.region_to_gene[assignment.edge_region_index] >= 0
    if not np.any(gene_edge):
        return np.zeros(0, dtype=np.int64)
    return np.unique(assignment.edge_col_index[gene_edge]).astype(np.int64)


def assignment_col_subset(assignment: CascadeAssignment, cols: np.ndarray) -> CascadeAssignment:
    """Restrict assignment edges and direct columns to a column index subset."""
    allowed = set(int(c) for c in np.asarray(cols, dtype=np.int64).tolist())
    if not allowed:
        raise ValueError("cols must be non-empty")
    if assignment.edge_col_index.size:
        edge_mask = np.asarray(
            [int(c) in allowed for c in assignment.edge_col_index.tolist()],
            dtype=bool,
        )
    else:
        edge_mask = np.zeros(0, dtype=bool)
    direct_mask = np.asarray(
        [int(c) in allowed for c in assignment.direct_col_index.tolist()],
        dtype=bool,
    )
    if edge_mask.size:
        edge_col = assignment.edge_col_index[edge_mask]
        edge_reg = assignment.edge_region_index[edge_mask]
        used_regions = np.unique(edge_reg.astype(np.int64, copy=False))
    else:
        edge_col = np.zeros(0, dtype=np.int64)
        edge_reg = np.zeros(0, dtype=np.int64)
        used_regions = np.zeros(0, dtype=np.int64)
    old_to_new = {int(old): i for i, old in enumerate(used_regions.tolist())}
    new_region_ids = [assignment.region_ids[int(i)] for i in used_regions.tolist()]
    new_region_type_id = assignment.region_type_id[used_regions] if used_regions.size else np.zeros(0, dtype=np.int64)
    new_region_to_gene = assignment.region_to_gene[used_regions] if used_regions.size else np.zeros(0, dtype=np.int64)
    new_orphan = new_region_to_gene < 0 if new_region_to_gene.size else np.zeros(0, dtype=bool)
    new_allocated = [assignment.allocated_gene_id[int(i)] for i in used_regions.tolist()]
    new_edge_reg = (
        np.asarray([old_to_new[int(r)] for r in edge_reg], dtype=np.int64) if edge_reg.size else edge_reg
    )
    new_direct = assignment.direct_col_index[direct_mask]
    return CascadeAssignment(
        gene_ids=list(assignment.gene_ids),
        region_ids=new_region_ids,
        region_type_id=new_region_type_id,
        region_to_gene=new_region_to_gene,
        orphan_region_mask=new_orphan,
        edge_col_index=edge_col.astype(np.int64, copy=False),
        edge_region_index=new_edge_reg,
        direct_col_index=new_direct.astype(np.int64, copy=False),
        region_types=assignment.region_types,
        n_study_loci=assignment.n_study_loci,
        allocated_gene_id=new_allocated,
    )


def assignment_gene_linked_only(assignment: CascadeAssignment) -> CascadeAssignment:
    """Restrict assignment to gene-linked typed edges; drop direct and orphan paths."""
    if assignment.edge_col_index.size == 0:
        return CascadeAssignment(
            gene_ids=list(assignment.gene_ids),
            region_ids=[],
            region_type_id=np.zeros(0, dtype=np.int64),
            region_to_gene=np.zeros(0, dtype=np.int64),
            orphan_region_mask=np.zeros(0, dtype=bool),
            edge_col_index=np.zeros(0, dtype=np.int64),
            edge_region_index=np.zeros(0, dtype=np.int64),
            direct_col_index=np.zeros(0, dtype=np.int64),
            region_types=assignment.region_types,
            n_study_loci=assignment.n_study_loci,
            allocated_gene_id=[],
        )
    gene_edge_mask = assignment.region_to_gene[assignment.edge_region_index] >= 0
    if not np.any(gene_edge_mask):
        raise ValueError("no gene-linked edges in assignment")
    edge_col = assignment.edge_col_index[gene_edge_mask]
    edge_reg = assignment.edge_region_index[gene_edge_mask]
    used_regions = np.unique(edge_reg.astype(np.int64, copy=False))
    old_to_new = {int(old): i for i, old in enumerate(used_regions.tolist())}
    new_region_ids = [assignment.region_ids[int(i)] for i in used_regions.tolist()]
    new_region_type_id = assignment.region_type_id[used_regions]
    new_region_to_gene = assignment.region_to_gene[used_regions]
    new_orphan_mask = new_region_to_gene < 0
    new_allocated = [assignment.allocated_gene_id[int(i)] for i in used_regions.tolist()]
    new_edge_region = np.asarray([old_to_new[int(r)] for r in edge_reg], dtype=np.int64)
    return CascadeAssignment(
        gene_ids=list(assignment.gene_ids),
        region_ids=new_region_ids,
        region_type_id=new_region_type_id,
        region_to_gene=new_region_to_gene,
        orphan_region_mask=new_orphan_mask,
        edge_col_index=edge_col.astype(np.int64, copy=False),
        edge_region_index=new_edge_region,
        direct_col_index=np.zeros(0, dtype=np.int64),
        region_types=assignment.region_types,
        n_study_loci=assignment.n_study_loci,
        allocated_gene_id=new_allocated,
    )


def _region_has_explicit_gene(gid: object) -> bool:
    if gid is None or (isinstance(gid, float) and np.isnan(gid)):
        return False
    return str(gid) not in ("", ".", "None", "nan")


def _midpoint_from_region(rec: object) -> tuple[str, int]:
    chrom = str(getattr(rec, "chromosome", ""))
    start = int(getattr(rec, "start", 0) or 0)
    end = int(getattr(rec, "end", start) or start)
    return chrom, (start + end) // 2


def _nearest_gene_distance(
    chromosome: str,
    midpoint: int,
    genes: pd.DataFrame,
) -> tuple[str | None, int | None]:
    """Return nearest gene_id on chromosome and distance in bp (0 if inside interval)."""
    if genes.empty or "chromosome" not in genes.columns:
        return None, None
    chrom = str(chromosome)
    g = genes.loc[genes["chromosome"].astype(str) == chrom]
    if g.empty:
        return None, None
    starts = g["start"].to_numpy(dtype=np.int64, copy=False)
    ends = g["end"].to_numpy(dtype=np.int64, copy=False)
    gene_ids = g["gene_id"].astype(str).to_numpy(copy=False)
    dist = np.where(
        (midpoint >= starts) & (midpoint <= ends),
        0,
        np.minimum(np.abs(midpoint - starts), np.abs(midpoint - ends)),
    )
    best = int(np.argmin(dist))
    return str(gene_ids[best]), int(dist[best])


def nearest_gene_on_chromosome(
    chromosome: str,
    midpoint: int,
    genes: pd.DataFrame,
) -> str | None:
    """Return nearest gene_id on ``chromosome`` by distance to gene interval; else None."""
    gid, _ = _nearest_gene_distance(chromosome, midpoint, genes)
    return gid


def _gene_tables_by_chrom(genes: pd.DataFrame) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """chromosome → (starts, ends, gene_ids) for fast nearest-gene."""
    out: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    if genes.empty or "chromosome" not in genes.columns:
        return out
    for chrom, g in genes.groupby(genes["chromosome"].astype(str), sort=False):
        out[str(chrom)] = (
            g["start"].to_numpy(dtype=np.int64, copy=False),
            g["end"].to_numpy(dtype=np.int64, copy=False),
            g["gene_id"].astype(str).to_numpy(copy=False),
        )
    return out


def _nearest_from_chrom_table(
    midpoint: int,
    starts: np.ndarray,
    ends: np.ndarray,
    gene_ids: np.ndarray,
) -> str | None:
    if starts.size == 0:
        return None
    dist = np.where(
        (midpoint >= starts) & (midpoint <= ends),
        0,
        np.minimum(np.abs(midpoint - starts), np.abs(midpoint - ends)),
    )
    return str(gene_ids[int(np.argmin(dist))])


def allocate_rbs_genes(
    regions: pd.DataFrame,
    genes: pd.DataFrame,
    *,
    policy: GeneAllocationPolicy = "legacy_nearest",
    max_nearest_gene_bp: int | None = None,
) -> pd.Series:
    """Map each region row to allocated gene_id per ``policy``.

    Regions already carrying a non-null ``gene_id`` keep it. Under ``legacy_nearest``,
    null-gene typed RBS regions receive same-chromosome nearest gene. Under
    ``explicit_only``, null-gene regions stay unallocated (orphan). Under
    ``bounded_nearest``, nearest gene is used only when distance ≤ ``max_nearest_gene_bp``.
    """
    if policy == "bounded_nearest" and max_nearest_gene_bp is None:
        raise ValueError("bounded_nearest requires max_nearest_gene_bp")
    chrom_tables = _gene_tables_by_chrom(genes)
    out: list[str | None] = []
    for rec in regions.itertuples(index=False):
        gid = getattr(rec, "gene_id", None)
        if _region_has_explicit_gene(gid):
            out.append(str(gid))
            continue
        if policy == "explicit_only":
            out.append(None)
            continue
        chrom, mid = _midpoint_from_region(rec)
        if policy == "bounded_nearest":
            nearest, dist = _nearest_gene_distance(chrom, mid, genes)
            if nearest is None or dist is None or dist > int(max_nearest_gene_bp or 0):
                out.append(None)
            else:
                out.append(nearest)
            continue
        table = chrom_tables.get(chrom)
        if table is None:
            out.append(None)
            continue
        out.append(_nearest_from_chrom_table(mid, table[0], table[1], table[2]))
    return pd.Series(out, index=regions.index, dtype=object)


def build_cascade_assignment(
    *,
    locus_index: pd.DataFrame,
    locus_region_edges: pd.DataFrame,
    regions: pd.DataFrame,
    genes: pd.DataFrame,
    max_loci: int | None = None,
    gene_allocation: GeneAllocationPolicy = "legacy_nearest",
    max_nearest_gene_bp: int | None = None,
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

    # Restrict edges + nearest-gene to regions that touch the study locus prefix.
    # Graph edges may store locus_id as float; normalize so "10.0" == "10".
    def _norm_ids(series: pd.Series) -> pd.Series:
        return series.map(
            lambda x: str(int(float(x))) if x is not None and not (isinstance(x, float) and np.isnan(x)) else ""
        )

    lr = locus_region_edges.loc[:, ["locus_id", "region_id"]].copy()
    lr["region_id"] = lr["region_id"].astype(str)
    lr["locus_id"] = _norm_ids(lr["locus_id"])
    study["locus_id"] = _norm_ids(study["locus_id"])
    study_locus_ids = set(study["locus_id"].tolist()) - {""}
    lr = lr.loc[lr["locus_id"].isin(study_locus_ids)].copy()
    touch_ids = set(lr["region_id"].tolist())
    selected = selected.loc[selected["region_id"].astype(str).isin(touch_ids)].copy()
    if selected.empty:
        # No typed edges in prefix → all direct.
        all_cols = study["col_index"].to_numpy(dtype=np.int64, copy=False)
        direct_col_index = all_cols[(all_cols >= 0) & (all_cols < n_study_loci)].astype(np.int64)
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

    selected["allocated_gene_id"] = allocate_rbs_genes(
        selected,
        genes,
        policy=gene_allocation,
        max_nearest_gene_bp=max_nearest_gene_bp,
    )

    region_view = selected.loc[
        :,
        [c for c in ("region_id", "region_type", "allocated_gene_id", "region_system") if c in selected.columns],
    ].copy()
    region_view["region_id"] = region_view["region_id"].astype(str)
    region_view["region_type"] = region_view["region_type"].astype(str)

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
