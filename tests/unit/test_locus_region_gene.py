"""Unit tests for hierarchical locus→region→gene index and packing."""

from __future__ import annotations

import pandas as pd
import torch

from mbs.models import HierarchicalDeepSet
from mbs.training.hier_dataset import (
    make_synthetic_hier_overfit_bundle,
    pack_hier_records_to_batch,
)
from mbs.training.locus_region_gene import (
    REGION_TYPE_TO_ID,
    UNASSIGNED_GENE_ID,
    UNASSIGNED_REGION_TYPE,
    build_locus_region_gene_index,
)
from mbs.training.multitask import MultitaskHeads, masked_multitask_loss


def test_build_locus_region_gene_index_retains_unassigned_singletons() -> None:
    locus_index = pd.DataFrame(
        {
            "col_index": [0, 1, 2],
            "locus_id": [10, 20, 30],
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
    assert index.n_unassigned_regions == 1
    assert UNASSIGNED_GENE_ID in index.gene_ids
    assert index.gene_ids[-1] == UNASSIGNED_GENE_ID
    orphan_region = f"unassigned:{30}"
    assert orphan_region in index.region_ids
    orphan_idx = index.region_ids.index(orphan_region)
    assert int(index.region_type_id[orphan_idx]) == REGION_TYPE_TO_ID[UNASSIGNED_REGION_TYPE]
    assert int(index.region_to_gene[orphan_idx]) == index.unassigned_gene_index
    # Typed roles preserved (not collapsed).
    assert set(index.region_type_id.tolist()) >= {
        REGION_TYPE_TO_ID["promoter_core"],
        REGION_TYPE_TO_ID["gene_body"],
        REGION_TYPE_TO_ID[UNASSIGNED_REGION_TYPE],
    }


def test_hier_pack_and_ablation_zeroes_disallowed_region_features() -> None:
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
    )
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
    )
    mbs = out["mbs"].view(len(records), batch.n_genes)
    present = out["present"].view(len(records), batch.n_genes)
    heads = MultitaskHeads(batch.n_genes, n_tissue_classes=3, sex_enabled=True)
    result = masked_multitask_loss(
        mbs=mbs,
        present=present,
        heads=heads,
        batch=batch,
    )
    assert torch.isfinite(result.loss)
    # Unassigned gene may be absent when only gene_body edges keep signal.
    if locus_region.unassigned_gene_index is not None:
        u = locus_region.unassigned_gene_index
        assert bool(present[:, u].eq(False).all())


def test_synthetic_hier_overfit_forward_shapes() -> None:
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
        n_region_types=6,
        dropout=0.0,
    )
    out = model(
        cpg_features=batch.cpg_features,
        cpg_to_region=batch.cpg_to_region,
        region_type=batch.region_type,
        region_to_gene=batch.region_to_gene,
        n_regions=3 * batch.n_regions,
        n_gene_instances=3 * batch.n_genes,
    )
    assert out["mbs"].shape == (3 * batch.n_genes,)
    assert out["region_present"].shape == (3 * batch.n_regions,)
