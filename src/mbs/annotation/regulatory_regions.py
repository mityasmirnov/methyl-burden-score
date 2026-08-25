"""Non-gene CGI/shore regions for graph-v2 RBS (ADR 0006)."""

from __future__ import annotations

from collections.abc import Collection

import duckdb
import pandas as pd

RBS_CONTEXTS = frozenset({"island", "north_shore", "south_shore"})

_EMPTY_REGIONS = pd.DataFrame(
    columns=[
        "region_id",
        "gene_id",
        "region_type",
        "chromosome",
        "start",
        "end",
        "strand",
        "source_version",
        "region_system",
    ]
)
_EMPTY_EDGES = pd.DataFrame(
    columns=["locus_id", "region_id", "edge_weight", "evidence_type", "primary_gene_role"]
)


def build_rbs_regions(
    loci: pd.DataFrame,
    *,
    gene_assigned_locus_ids: Collection[object],
    islands: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """RBS regions from CGI/shore context on loci not in gene five-role edges.

    When ``islands`` is provided (UCSC-style ``chromosome,start,end``), each
    island interval gets distinct island/shore regions keyed to the nearest
    island. Without islands, fall back to one region per chromosome × context
    (fixture / no-CGI path).
    """
    assigned = set(gene_assigned_locus_ids)
    if loci.empty or "cpg_context" not in loci.columns:
        return _EMPTY_REGIONS.copy(), _EMPTY_EDGES.copy()
    cand = loci.loc[
        (~loci["locus_id"].isin(list(assigned))) & loci["cpg_context"].isin(list(RBS_CONTEXTS))
    ].copy()
    if cand.empty:
        return _EMPTY_REGIONS.copy(), _EMPTY_EDGES.copy()
    if islands is not None and not islands.empty:
        return _build_per_island_rbs(cand, islands)
    return _build_chrom_context_rbs(cand)


def _build_chrom_context_rbs(cand: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """One RBS region per chromosome × CGI context (fixture fallback)."""
    rows: list[dict[str, object]] = []
    edges: list[dict[str, object]] = []
    for key, group in cand.groupby(["chromosome", "cpg_context"], sort=True):
        chrom, ctx = key  # type: ignore[misc]
        ordered = group.sort_values("position")
        start = int(ordered["position"].min())  # type: ignore[arg-type]
        end = int(ordered["position"].max())  # type: ignore[arg-type]
        region_id = f"rbs:{ctx}:{chrom}:{start}-{end}"
        rows.append(
            {
                "region_id": region_id,
                "gene_id": None,
                "region_type": f"cgi_{ctx}",
                "chromosome": chrom,
                "start": start,
                "end": end,
                "strand": ".",
                "source_version": "UCSC_cgi_context",
                "region_system": "rbs",
            }
        )
        edges.extend(
            {
                "locus_id": lid,
                "region_id": region_id,
                "edge_weight": 1.0,
                "evidence_type": "cgi_unassigned",
                "primary_gene_role": False,
            }
            for lid in ordered["locus_id"].tolist()
        )
    return pd.DataFrame(rows), pd.DataFrame(edges)


def _build_per_island_rbs(
    cand: pd.DataFrame,
    islands: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """One island + shore region per UCSC island, nearest-island assignment."""
    required = {"chromosome", "start", "end"}
    missing = required - set(islands.columns)
    if missing:
        raise ValueError(f"islands missing columns: {sorted(missing)}")
    con = duckdb.connect(database=":memory:")
    con.register(
        "cand",
        cand[["locus_id", "chromosome", "position", "cpg_context"]],
    )
    con.register("islands", islands[["chromosome", "start", "end"]])
    assigned = con.execute(
        """
        WITH dist AS (
          SELECT
            c.locus_id,
            c.chromosome,
            c.position,
            c.cpg_context,
            i.start AS island_start,
            i.end AS island_end,
            CASE
              WHEN c.position BETWEEN i.start AND i.end THEN 0
              WHEN c.position < i.start THEN i.start - c.position
              ELSE c.position - i.end
            END AS abs_dist
          FROM cand c
          JOIN islands i
            ON c.chromosome = i.chromosome
        ),
        best AS (
          SELECT
            *,
            ROW_NUMBER() OVER (
              PARTITION BY locus_id
              ORDER BY abs_dist ASC, island_start ASC, island_end ASC
            ) AS rn
          FROM dist
        )
        SELECT
          locus_id,
          chromosome,
          position,
          cpg_context,
          island_start,
          island_end
        FROM best
        WHERE rn = 1
        """
    ).fetchdf()
    con.close()
    if assigned.empty:
        return _EMPTY_REGIONS.copy(), _EMPTY_EDGES.copy()

    rows: list[dict[str, object]] = []
    edges: list[dict[str, object]] = []
    seen_regions: set[str] = set()
    for key, group in assigned.groupby(
        ["chromosome", "cpg_context", "island_start", "island_end"],
        sort=True,
    ):
        chrom, ctx, i_start, i_end = key  # type: ignore[misc]
        island_start = int(i_start)  # type: ignore[arg-type]
        island_end = int(i_end)  # type: ignore[arg-type]
        region_id = f"rbs:{ctx}:{chrom}:{island_start}-{island_end}"
        if region_id not in seen_regions:
            seen_regions.add(region_id)
            if ctx == "island":
                start, end = island_start, island_end
            else:
                start = int(group["position"].min())  # type: ignore[arg-type]
                end = int(group["position"].max())  # type: ignore[arg-type]
            rows.append(
                {
                    "region_id": region_id,
                    "gene_id": None,
                    "region_type": f"cgi_{ctx}",
                    "chromosome": chrom,
                    "start": start,
                    "end": end,
                    "strand": ".",
                    "source_version": "UCSC_cgi_per_island_shore",
                    "region_system": "rbs",
                }
            )
        edges.extend(
            {
                "locus_id": lid,
                "region_id": region_id,
                "edge_weight": 1.0,
                "evidence_type": "cgi_unassigned",
                "primary_gene_role": False,
            }
            for lid in group["locus_id"].tolist()
        )
    return pd.DataFrame(rows), pd.DataFrame(edges)
