"""Non-gene CGI/shore regions for graph-v2 RBS (ADR 0006)."""

from __future__ import annotations

from collections.abc import Collection

import pandas as pd

RBS_CONTEXTS = frozenset({"island", "north_shore", "south_shore"})


def build_rbs_regions(
    loci: pd.DataFrame,
    *,
    gene_assigned_locus_ids: Collection[object],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """RBS regions from CGI/shore context on loci not in gene five-role edges."""
    assigned = set(gene_assigned_locus_ids)
    if loci.empty or "cpg_context" not in loci.columns:
        empty_r = pd.DataFrame(
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
        empty_e = pd.DataFrame(
            columns=["locus_id", "region_id", "edge_weight", "evidence_type", "primary_gene_role"]
        )
        return empty_r, empty_e
    cand = loci.loc[
        (~loci["locus_id"].isin(assigned)) & loci["cpg_context"].isin(RBS_CONTEXTS)
    ].copy()
    if cand.empty:
        empty_r = pd.DataFrame(
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
        empty_e = pd.DataFrame(
            columns=["locus_id", "region_id", "edge_weight", "evidence_type", "primary_gene_role"]
        )
        return empty_r, empty_e
    rows: list[dict[str, object]] = []
    edges: list[dict[str, object]] = []
    for (chrom, ctx), group in cand.groupby(["chromosome", "cpg_context"], sort=True):
        ordered = group.sort_values("position")
        start = int(ordered["position"].min())
        end = int(ordered["position"].max())
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
        for lid in ordered["locus_id"]:
            edges.append(  # noqa: PERF401
                {
                    "locus_id": lid,
                    "region_id": region_id,
                    "edge_weight": 1.0,
                    "evidence_type": "cgi_unassigned",
                    "primary_gene_role": False,
                }
            )
    return pd.DataFrame(rows), pd.DataFrame(edges)
