"""Unit checks for annotation coverage aggregation."""

from __future__ import annotations

import pandas as pd

from mbs.annotation.coverage import compute_annotation_coverage


def test_compute_annotation_coverage_tiny() -> None:
    probes = pd.DataFrame(
        {
            "probe_id": ["a", "b", "c", "d"],
            "platform_id": ["HM450", "HM450", "HM450", "EPIC"],
            "mapping_status": ["mapped", "mapped", "unmapped", "mapped"],
        }
    )
    probe_locus_edges = pd.DataFrame(
        {
            "probe_id": ["a", "b", "d"],
            "platform_id": ["HM450", "HM450", "EPIC"],
            "locus_id": [1, 2, 1],
        }
    )
    loci = pd.DataFrame({"locus_id": [1, 2, 3]})
    regions = pd.DataFrame(
        {
            "region_id": ["r1", "r2"],
            "region_type": ["promoter_core", "gene_body"],
            "gene_id": ["g1", "g1"],
        }
    )
    locus_region_edges = pd.DataFrame(
        {
            "locus_id": [1],
            "region_id": ["r1"],
            "edge_weight": [1.0],
            "evidence_type": ["overlap"],
            "primary_gene_role": ["promoter_core"],
        }
    )
    out = compute_annotation_coverage(
        probes=probes,
        probe_locus_edges=probe_locus_edges,
        loci=loci,
        locus_region_edges=locus_region_edges,
        regions=regions,
        graph_id="test-graph",
    )
    assert out["locus_level"]["n_loci"] == 3
    assert out["locus_level"]["n_assigned_loci"] == 1
    assert out["locus_level"]["n_unassigned_loci"] == 2
    hm450 = next(p for p in out["probe_level"]["platforms"] if p["platform_id"] == "HM450")
    assert hm450["n_probes"] == 3
    assert hm450["n_mapped"] == 2
    assert hm450["n_unmapped"] == 1
    assert hm450["n_mapped_assigned"] == 1  # probe a → locus 1
    assert hm450["n_mapped_unassigned"] == 1  # probe b → locus 2
    assert hm450["probes_by_role"]["promoter_core"] == 1
    epic = next(p for p in out["probe_level"]["platforms"] if p["platform_id"] == "EPIC")
    assert epic["n_mapped_assigned"] == 1
    assert epic["n_unmapped"] == 0
