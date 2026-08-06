"""Build the shared GRCh38 locus registry from platform probes."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from mbs.annotation.cgi_context import annotate_cpg_context
from mbs.annotation.export_infinium import load_infinium_probes


def _canonical_key(chromosome: str, position: int) -> str:
    return f"GRCh38:{chromosome}:{position}"


def build_locus_registry(
    probes: pd.DataFrame,
    cgi_path: Path | None = None,
    *,
    shore_bp: int = 2000,
    shelf_bp: int = 4000,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Collapse probes to canonical loci and emit probe/locus edge tables.

    Returns ``(loci, probes_out, probe_locus_edges)``.
    """
    required = {
        "probe_id",
        "platform_id",
        "chromosome",
        "position",
        "mapping_status",
        "M_mapping",
        "M_nonuniq",
        "M_general",
    }
    missing = required - set(probes.columns)
    if missing:
        raise KeyError(f"probe table missing columns: {sorted(missing)}")

    mapped = probes[probes["mapping_status"] == "mapped"].copy()
    locus_keys = mapped.loc[:, ["chromosome", "position"]].drop_duplicates()
    locus_keys = locus_keys.sort_values(
        by=["chromosome", "position"],
        kind="mergesort",
    ).reset_index(drop=True)
    locus_keys["locus_id"] = (locus_keys.index + 1).astype("uint64")
    locus_keys["genome_build"] = "GRCh38"
    locus_keys["canonical_key"] = [
        _canonical_key(str(chrom), int(pos))
        for chrom, pos in zip(locus_keys["chromosome"], locus_keys["position"], strict=True)
    ]
    locus_keys["mapping_status"] = "mapped"

    if cgi_path is not None:
        locus_keys = annotate_cpg_context(
            locus_keys,
            cgi_path,
            shore_bp=shore_bp,
            shelf_bp=shelf_bp,
        )
    else:
        locus_keys["cpg_context"] = "open_sea"

    loci = locus_keys.loc[
        :,
        [
            "locus_id",
            "genome_build",
            "chromosome",
            "position",
            "canonical_key",
            "mapping_status",
            "cpg_context",
        ],
    ].copy()

    probe_cols = [
        "probe_id",
        "platform_id",
        "probe_design",
        "core_probe_id",
        "M_mapping",
        "M_nonuniq",
        "M_general",
        "mapQ",
        "strand",
        "mapping_status",
    ]
    present_cols = [c for c in probe_cols if c in probes.columns]
    probes_out = probes.loc[:, present_cols].copy()

    edges = mapped.merge(
        loci.loc[:, ["locus_id", "chromosome", "position"]],
        on=["chromosome", "position"],
        how="left",
        validate="many_to_one",
    )
    edges["mapping_source"] = "infinium_annotation_hg38"
    edges["mapping_confidence"] = edges["mapQ"].astype("Float64") / 60.0
    edges["is_primary"] = True
    probe_locus_edges = edges.loc[
        :,
        [
            "probe_id",
            "platform_id",
            "locus_id",
            "mapping_source",
            "mapping_confidence",
            "is_primary",
        ],
    ].copy()
    probe_locus_edges["locus_id"] = probe_locus_edges["locus_id"].astype("uint64")
    return loci, probes_out, probe_locus_edges


def build_locus_registry_from_infinium(
    infinium_root: Path,
    platforms: tuple[str, ...] | list[str],
    cgi_path: Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load Infinium tables and build the registry.

    Returns ``(loci, probes, edges, raw_probes)``.
    """
    raw_probes = load_infinium_probes(infinium_root, platforms)
    loci, probes, edges = build_locus_registry(raw_probes, cgi_path=cgi_path)
    return loci, probes, edges, raw_probes
