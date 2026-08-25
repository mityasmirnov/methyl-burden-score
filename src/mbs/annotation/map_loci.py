"""Map canonical loci to gene regions with Stage 0 role precedence."""

from __future__ import annotations

from pathlib import Path

import duckdb
import pandas as pd

from mbs.annotation.gencode_regions import ROLE_PRECEDENCE

_PRECEDENCE_CASE = " ".join(
    f"WHEN region_type = '{role}' THEN {rank}" for rank, role in enumerate(ROLE_PRECEDENCE)
)


def map_loci_to_regions(
    loci: pd.DataFrame,
    regions: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return ``(locus_region_edges, region_gene_edges)``.

    For each ``(locus_id, gene_id)`` pair, keep only the highest-precedence
    overlapping region type. A locus may still map to multiple genes.
    """
    if loci.empty or regions.empty:
        empty_lr = pd.DataFrame(
            columns=[
                "locus_id",
                "region_id",
                "edge_weight",
                "evidence_type",
                "primary_gene_role",
            ]
        )
        empty_rg = pd.DataFrame(columns=["region_id", "gene_id", "edge_weight"])
        return empty_lr, empty_rg

    con = duckdb.connect(database=":memory:")
    con.register("loci", loci[["locus_id", "chromosome", "position"]])
    con.register(
        "regions",
        regions[["region_id", "gene_id", "region_type", "chromosome", "start", "end"]],
    )
    edges = con.execute(
        f"""
        WITH hits AS (
          SELECT
            l.locus_id,
            r.region_id,
            r.gene_id,
            r.region_type,
            CASE {_PRECEDENCE_CASE} ELSE 999 END AS precedence_rank
          FROM loci l
          JOIN regions r
            ON l.chromosome = r.chromosome
           AND l.position BETWEEN r.start AND r.end
        ),
        ranked AS (
          SELECT
            *,
            ROW_NUMBER() OVER (
              PARTITION BY locus_id, gene_id
              ORDER BY precedence_rank ASC, region_id ASC
            ) AS rn
          FROM hits
        )
        SELECT
          locus_id,
          region_id,
          1.0 AS edge_weight,
          'interval_overlap' AS evidence_type,
          TRUE AS primary_gene_role
        FROM ranked
        WHERE rn = 1
        ORDER BY locus_id, region_id
        """
    ).fetchdf()
    con.close()

    edges["locus_id"] = edges["locus_id"].astype("uint64")
    edges["edge_weight"] = edges["edge_weight"].astype("float64")
    edges["primary_gene_role"] = True

    region_gene = (
        regions[["region_id", "gene_id"]]
        .drop_duplicates()
        .assign(edge_weight=1.0)
        .sort_values(["gene_id", "region_id"], kind="mergesort")
        .reset_index(drop=True)
    )
    return edges, region_gene


def write_regions_bed(regions: pd.DataFrame, path: Path) -> None:
    """Write BED6+ with gene_id and region_type columns (null gene_id → ``.``)."""
    gene_col = regions["gene_id"].where(regions["gene_id"].notna(), ".")
    bed = pd.DataFrame(
        {
            "chrom": regions["chromosome"],
            "start": regions["start"].astype("int64") - 1,  # 0-based BED
            "end": regions["end"].astype("int64"),
            "region_id": regions["region_id"],
            "score": 0,
            "strand": regions["strand"],
            "gene_id": gene_col,
            "region_type": regions["region_type"],
        }
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    bed.to_csv(path, sep="\t", header=False, index=False)
