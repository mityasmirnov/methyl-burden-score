"""Per-platform probe→locus→region coverage joins for Stage 0 annotation."""

from __future__ import annotations

from typing import Any

import pandas as pd

ROLE_ORDER = (
    "promoter_core",
    "promoter_proximal",
    "five_prime",
    "gene_body",
    "three_prime",
)


def compute_annotation_coverage(
    *,
    probes: pd.DataFrame,
    probe_locus_edges: pd.DataFrame,
    loci: pd.DataFrame,
    locus_region_edges: pd.DataFrame,
    regions: pd.DataFrame,
    graph_id: str,
) -> dict[str, Any]:
    """Join Illumina probes to five-role regions; return JSON-serializable stats.

    Distinguishes:
    - unmapped probe: ``mapping_status != mapped`` (no GRCh38 cytosine)
    - unassigned locus/probe: mapped but no ``locus_region_edges`` row
    """
    required_probe = {"probe_id", "platform_id", "mapping_status"}
    missing = required_probe - set(probes.columns)
    if missing:
        raise KeyError(f"probes missing columns: {sorted(missing)}")

    assigned_locus_ids = set(locus_region_edges["locus_id"].astype("uint64").tolist())
    loci_ids = set(loci["locus_id"].astype("uint64").tolist())
    n_loci = len(loci)
    n_unassigned_loci = len(loci_ids - assigned_locus_ids)
    n_assigned_loci = n_loci - n_unassigned_loci

    lr = locus_region_edges.merge(
        regions[["region_id", "region_type"]],
        on="region_id",
        how="left",
        validate="many_to_one",
    )
    locus_by_role = {
        role: int(lr.loc[lr["region_type"] == role, "locus_id"].nunique()) for role in ROLE_ORDER
    }
    edge_by_role = {role: int((lr["region_type"] == role).sum()) for role in ROLE_ORDER}

    multi_gene = 0
    if not lr.empty:
        genes_per_locus = (
            lr.merge(
                regions[["region_id", "gene_id"]],
                on="region_id",
                how="left",
            )
            .groupby("locus_id")["gene_id"]
            .nunique()
        )
        multi_gene = int((genes_per_locus > 1).sum())

    platform_rows: list[dict[str, Any]] = []
    for platform_id, psub in probes.groupby("platform_id", sort=True):
        n_probes_plat = len(psub)
        n_mapped = int((psub["mapping_status"] == "mapped").sum())
        n_unmapped = n_probes_plat - n_mapped
        pedges = probe_locus_edges.loc[
            probe_locus_edges["platform_id"] == platform_id,
            ["probe_id", "locus_id"],
        ].copy()
        pedges["locus_id"] = pedges["locus_id"].astype("uint64")
        pedges["assigned"] = pedges["locus_id"].isin(assigned_locus_ids)
        n_mapped_assigned = int(pedges["assigned"].sum())
        n_mapped_unassigned = int((~pedges["assigned"]).sum())
        joined = pedges.merge(
            lr[["locus_id", "region_type"]],
            on="locus_id",
            how="left",
        )
        probe_by_role = {
            role: int(joined.loc[joined["region_type"] == role, "probe_id"].nunique())
            for role in ROLE_ORDER
        }
        platform_rows.append(
            {
                "platform_id": str(platform_id),
                "n_probes": n_probes_plat,
                "n_mapped": n_mapped,
                "n_unmapped": n_unmapped,
                "pct_unmapped": (
                    round(100.0 * n_unmapped / n_probes_plat, 4) if n_probes_plat else 0.0
                ),
                "n_mapped_assigned": n_mapped_assigned,
                "n_mapped_unassigned": n_mapped_unassigned,
                "pct_mapped_assigned": (
                    round(100.0 * n_mapped_assigned / n_mapped, 4) if n_mapped else 0.0
                ),
                "pct_mapped_unassigned": (
                    round(100.0 * n_mapped_unassigned / n_mapped, 4) if n_mapped else 0.0
                ),
                "probes_by_role": probe_by_role,
            }
        )

    n_probes = len(probes)
    n_unmapped_probes = int((probes["mapping_status"] != "mapped").sum())
    return {
        "graph_id": graph_id,
        "locus_level": {
            "n_loci": n_loci,
            "n_assigned_loci": n_assigned_loci,
            "n_unassigned_loci": n_unassigned_loci,
            "pct_assigned": round(100.0 * n_assigned_loci / n_loci, 4) if n_loci else 0.0,
            "pct_unassigned": round(100.0 * n_unassigned_loci / n_loci, 4) if n_loci else 0.0,
            "n_loci_multi_gene": multi_gene,
            "loci_by_role": locus_by_role,
            "locus_region_edges_by_role": edge_by_role,
            "n_locus_region_edges": len(locus_region_edges),
        },
        "probe_level": {
            "n_probes": n_probes,
            "n_unmapped_probes": n_unmapped_probes,
            "pct_unmapped": (round(100.0 * n_unmapped_probes / n_probes, 4) if n_probes else 0.0),
            "n_probe_locus_edges": len(probe_locus_edges),
            "platforms": platform_rows,
        },
        "definitions": {
            "unmapped_probe": (
                "Illumina probe with no GRCh38 cytosine coordinate; excluded from matrices"
            ),
            "unassigned_locus": (
                "Mapped locus with no five-role GENCODE region edge; "
                "hier train uses singleton unassigned"
            ),
            "atlas_probe_annotations": (
                "EWAS Atlas probe TSV is a separate unused layer; not counted here"
            ),
        },
    }
