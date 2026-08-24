"""Adaptive CpG-count tiles for leftover mapped loci (TBS)."""

from __future__ import annotations

from collections.abc import Collection

import pandas as pd

DEFAULT_TILE_N = 50


def build_tiles(
    loci: pd.DataFrame,
    *,
    remaining_locus_ids: Collection[object],
    target_n_cpgs: int = DEFAULT_TILE_N,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Greedy chrom-sorted tiles of ``target_n_cpgs`` loci."""
    if target_n_cpgs < 1:
        raise ValueError("target_n_cpgs must be >= 1")
    remain = set(remaining_locus_ids)
    sub = loci.loc[loci["locus_id"].isin(remain)].copy()
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
    if sub.empty:
        return empty_r, empty_e
    sub = sub.sort_values(["chromosome", "position"])
    rows: list[dict[str, object]] = []
    edges: list[dict[str, object]] = []
    tile_i = 0
    for chrom, group in sub.groupby("chromosome", sort=True):
        chunk: list[tuple[object, int]] = []
        for rec in group.itertuples(index=False):
            chunk.append((rec.locus_id, int(rec.position)))
            if len(chunk) >= target_n_cpgs:
                _emit_tile(rows, edges, chrom, tile_i, chunk)
                tile_i += 1
                chunk = []
        if chunk:
            _emit_tile(rows, edges, chrom, tile_i, chunk)
            tile_i += 1
    return pd.DataFrame(rows), pd.DataFrame(edges)


def _emit_tile(
    rows: list[dict[str, object]],
    edges: list[dict[str, object]],
    chrom: object,
    tile_i: int,
    chunk: list[tuple[object, int]],
) -> None:
    start = chunk[0][1]
    end = chunk[-1][1]
    region_id = f"tbs:{chrom}:t{tile_i}"
    rows.append(
        {
            "region_id": region_id,
            "gene_id": None,
            "region_type": "cpg_tile",
            "chromosome": chrom,
            "start": start,
            "end": end,
            "strand": ".",
            "source_version": "adaptive_cpg_count",
            "region_system": "tbs",
        }
    )
    for lid, _pos in chunk:
        edges.append(
            {
                "locus_id": lid,
                "region_id": region_id,
                "edge_weight": 1.0,
                "evidence_type": "tile",
                "primary_gene_role": False,
            }
        )
