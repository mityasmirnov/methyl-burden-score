"""Unit tests for hierarchical locus→region→gene index and residual path."""

from __future__ import annotations

import pandas as pd
import torch

from mbs.annotation.build import attach_cgi_tile_systems
from mbs.batch import (
    ANNOTATION_STATUS_MAPPED,
    ANNOTATION_STATUS_UNMAPPED,
)
from mbs.models import HierarchicalDeepSet
from mbs.training.hier_dataset import (
    make_synthetic_hier_overfit_bundle,
    pack_hier_records_to_batch,
)
from mbs.training.locus_gene import build_locus_gene_index
from mbs.training.locus_region_gene import (
    REGION_TYPE_TO_ID,
    RESIDUAL_PANEL_ID,
    build_locus_region_gene_index,
)
from mbs.training.multitask import MultitaskHeads, masked_multitask_loss


def test_build_locus_region_gene_index_routes_orphans_to_residual() -> None:
    locus_index = pd.DataFrame(
        {
            "col_index": [0, 1, 2],
            "locus_id": [10, 20, 30],
            "canonical_key": ["chr1:10", "chr1:20", "chr1:30"],
        }
    )
    locus_region_edges = pd.DataFrame(
        {
            "locus_id": [10, 20],
            "region_id": ["G1:promoter_core", "G1:gene_body"],
        }
    )
    regions = pd.DataFrame(
        {
            "region_id": ["G1:promoter_core", "G1:gene_body"],
            "gene_id": ["G1", "G1"],
            "region_type": ["promoter_core", "gene_body"],
        }
    )
    index = build_locus_region_gene_index(
        locus_index=locus_index,
        locus_region_edges=locus_region_edges,
        regions=regions,
    )
    assert index.n_typed_edges == 2
    assert index.n_residual_cols == 1
    assert RESIDUAL_PANEL_ID not in index.gene_ids
    assert int(index.residual_col_index[0]) == 2
    assert int(index.column_annotation_status[0]) == ANNOTATION_STATUS_MAPPED
    assert int(index.column_annotation_status[2]) == ANNOTATION_STATUS_UNMAPPED
    assert set(index.region_type_id.tolist()) >= {
        REGION_TYPE_TO_ID["promoter_core"],
        REGION_TYPE_TO_ID["gene_body"],
    }


def test_hier_pack_mapped_only_drops_residual() -> None:
    bundle = make_synthetic_hier_overfit_bundle(n_samples=4, seed=3)
    locus_region = bundle["locus_region"]
    records = bundle["records"][:2]
    body_only = {REGION_TYPE_TO_ID["gene_body"]}
    batch = pack_hier_records_to_batch(
        records,
        locus_region=locus_region,
        age_values=[20.0, 30.0],
        age_enabled=[True, True],
        tissue_enabled=[True, True],
        sex_enabled=[True, False],
        sex_class_indices=[0, 1],
        allowed_region_type_ids=body_only,
        include_residual=False,
    )
    assert batch.residual_features.shape[0] == 0
    model = HierarchicalDeepSet(
        int(bundle["input_dim"]),
        n_region_types=len(bundle["region_types"]),
        dropout=0.0,
    )
    model.eval()
    out = model(
        cpg_features=batch.cpg_features,
        cpg_to_region=batch.cpg_to_region,
        region_type=batch.region_type,
        region_to_gene=batch.region_to_gene,
        n_regions=len(records) * batch.n_regions,
        n_gene_instances=len(records) * batch.n_genes,
        residual_features=batch.residual_features,
        residual_sample_index=batch.residual_sample_index,
        n_samples=len(records),
    )
    gene_mbs = out["mbs"].view(len(records), batch.n_genes)
    gene_present = out["present"].view(len(records), batch.n_genes)
    residual_mbs = out["residual_mbs"].view(len(records), 1)
    residual_present = out["residual_present"].view(len(records), 1)
    mbs = torch.cat([gene_mbs, residual_mbs], dim=1)
    present = torch.cat([gene_present, residual_present], dim=1)
    heads = MultitaskHeads(batch.n_genes + 1, n_tissue_classes=3, sex_enabled=True)
    result = masked_multitask_loss(
        mbs=mbs,
        present=present,
        heads=heads,
        batch=batch,
    )
    assert torch.isfinite(result.loss)
    assert bool(residual_present.eq(False).all())


def test_synthetic_hier_overfit_forward_shapes_with_residual() -> None:
    bundle = make_synthetic_hier_overfit_bundle(n_samples=3, seed=1)
    locus_region = bundle["locus_region"]
    batch = pack_hier_records_to_batch(
        bundle["records"],
        locus_region=locus_region,
        age_values=list(bundle["ages"]),
        age_enabled=[True] * 3,
        tissue_enabled=[True] * 3,
        sex_enabled=[True] * 3,
        sex_class_indices=[0, 1, 0],
    )
    model = HierarchicalDeepSet(
        int(bundle["input_dim"]),
        n_region_types=5,
        dropout=0.0,
    )
    out = model(
        cpg_features=batch.cpg_features,
        cpg_to_region=batch.cpg_to_region,
        region_type=batch.region_type,
        region_to_gene=batch.region_to_gene,
        n_regions=3 * batch.n_regions,
        n_gene_instances=3 * batch.n_genes,
        residual_features=batch.residual_features,
        residual_sample_index=batch.residual_sample_index,
        n_samples=3,
    )
    assert out["mbs"].shape == (3 * batch.n_genes,)
    assert out["region_present"].shape == (3 * batch.n_regions,)
    assert out["residual_mbs"].shape == (3,)
    assert out["residual_present"].shape == (3,)
    assert batch.annotation_status.numel() == batch.cpg_features.shape[0]


def test_multi_system_index_includes_rbs_tbs_gene_only_ignores() -> None:
    locus_index = pd.DataFrame(
        {
            "col_index": [0, 1, 2, 3],
            "locus_id": [10, 20, 30, 40],
            "canonical_key": ["chr1:10", "chr1:20", "chr1:30", "chr1:40"],
        }
    )
    regions = pd.DataFrame(
        {
            "region_id": ["G1:gene_body", "rbs:island:chr1:1-5", "tbs:chr1:t0"],
            "gene_id": ["G1", None, None],
            "region_type": ["gene_body", "cgi_island", "cpg_tile"],
            "region_system": ["gene", "rbs", "tbs"],
        }
    )
    edges = pd.DataFrame(
        {
            "locus_id": [10, 20, 30],
            "region_id": ["G1:gene_body", "rbs:island:chr1:1-5", "tbs:chr1:t0"],
        }
    )
    gene_only = build_locus_region_gene_index(
        locus_index=locus_index,
        locus_region_edges=edges,
        regions=regions,
        region_systems=("gene",),
    )
    assert gene_only.gene_ids == ["G1"]
    assert gene_only.n_typed_edges == 1
    assert set(gene_only.residual_col_index.tolist()) == {1, 2, 3}

    multi = build_locus_region_gene_index(
        locus_index=locus_index,
        locus_region_edges=edges,
        regions=regions,
        region_systems=("gene", "rbs", "tbs"),
    )
    assert set(multi.gene_ids) == {"G1", "rbs:island:chr1:1-5", "tbs:chr1:t0"}
    assert multi.n_typed_edges == 3
    assert multi.residual_col_index.tolist() == [3]
    assert "cgi_island" in multi.region_types
    assert "cpg_tile" in multi.region_types


def test_flat_arm_panels_disjoint() -> None:
    loci = pd.DataFrame(
        {
            "locus_id": [1, 2, 3, 4, 5, 6],
            "chromosome": ["chr1"] * 6,
            "position": [10, 20, 30, 40, 50, 60],
            "cpg_context": ["island", "island", "open_sea", "open_sea", "open_sea", "open_sea"],
            "mapping_status": ["mapped"] * 5 + ["unmapped"],
        }
    )
    regions = pd.DataFrame(
        {
            "region_id": ["g:body"],
            "gene_id": ["ENSG1"],
            "region_type": ["gene_body"],
            "chromosome": ["chr1"],
            "start": [1],
            "end": [15],
            "strand": ["+"],
            "source_version": ["gencode"],
            "region_system": ["gene"],
        }
    )
    edges = pd.DataFrame(
        {
            "locus_id": [1],
            "region_id": ["g:body"],
            "edge_weight": [1.0],
            "evidence_type": ["gene"],
            "primary_gene_role": [True],
        }
    )
    out_r, out_e = attach_cgi_tile_systems(loci, regions, edges, tile_target_n_cpgs=2)
    locus_index = pd.DataFrame({"col_index": list(range(6)), "locus_id": loci["locus_id"]})
    gene_idx = build_locus_gene_index(
        locus_index=locus_index,
        locus_region_edges=out_e,
        regions=out_r,
        region_systems=("gene",),
    )
    rbs_idx = build_locus_gene_index(
        locus_index=locus_index,
        locus_region_edges=out_e,
        regions=out_r,
        region_systems=("rbs",),
    )
    tbs_idx = build_locus_gene_index(
        locus_index=locus_index,
        locus_region_edges=out_e,
        regions=out_r,
        region_systems=("tbs",),
    )
    assert set(gene_idx.gene_ids).isdisjoint(rbs_idx.gene_ids)
    assert set(gene_idx.gene_ids).isdisjoint(tbs_idx.gene_ids)
    assert set(rbs_idx.gene_ids).isdisjoint(tbs_idx.gene_ids)
    gene_cols = set(gene_idx.edge_col_index.tolist())
    rbs_cols = set(rbs_idx.edge_col_index.tolist())
    tbs_cols = set(tbs_idx.edge_col_index.tolist())
    assert gene_cols.isdisjoint(rbs_cols)
    assert gene_cols.isdisjoint(tbs_cols)
    assert rbs_cols.isdisjoint(tbs_cols)
