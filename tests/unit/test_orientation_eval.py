"""Orientation contract v2 and flat_region feature path tests."""

from __future__ import annotations

import numpy as np
import torch

from mbs.scoring.orientation import (
    accumulate_signed_gene_mean_m,
    orient_mbs_array,
    signed_gene_mean_m,
)
from mbs.training.feature_schema import FLAT_REGION, m_column_index, observed_column_index
from mbs.training.flat_region_features import (
    apply_flat_region_feature_mode,
    flat_region_input_dim,
    gather_flat_region_features,
)
from mbs.training.multitask import MultitaskHeads


def test_flat_region_m_column_is_zero() -> None:
    assert m_column_index(feature_schema=FLAT_REGION, include_m_value=True) == 0
    assert observed_column_index(feature_schema=FLAT_REGION) == -1


def test_observed_mask_excludes_unobserved_from_signed_m() -> None:
    m = np.array([1.0, 99.0, 3.0], dtype=np.float64)
    genes = np.array([0, 0, 1], dtype=np.int64)
    obs = np.array([True, False, True], dtype=bool)
    out = signed_gene_mean_m(m, genes, n_genes=2, observed_mask=obs)
    np.testing.assert_allclose(out, [1.0, 3.0])


def test_accumulate_signed_gene_mean_m_observed_batches() -> None:
    out = accumulate_signed_gene_mean_m(
        n_genes=2,
        cpg_m_batches=[np.array([1.0, 5.0]), np.array([3.0])],
        cpg_to_gene_batches=[np.array([0, 0]), np.array([1])],
        observed_batches=[np.array([True, False]), np.array([True])],
    )
    np.testing.assert_allclose(out, [1.0, 3.0])


def test_orient_mbs_array_flips() -> None:
    mbs = np.array([[0.2, 0.8]], dtype=np.float32)
    flipped = orient_mbs_array(mbs, score_polarity="flipped")
    np.testing.assert_allclose(flipped, [[0.8, 0.2]])


def test_multitask_head_logit_invariance_under_orientation() -> None:
    """W·(x-0.5) equals (-W)·((1-x)-0.5) with unchanged bias (ADR 0008 pairing)."""
    n_genes = 4
    heads = MultitaskHeads(n_genes, n_tissue_classes=3, sex_enabled=False)
    with torch.no_grad():
        heads.tissue_head.gene_weight.normal_(0.0, 0.1)
        heads.tissue_head.bias.zero_()
    mbs = torch.rand(2, n_genes)
    present = torch.ones(2, n_genes, dtype=torch.bool)
    logits_raw = heads.forward_tissue(mbs, present)
    with torch.no_grad():
        w = heads.tissue_head.gene_weight.clone()
        heads.tissue_head.gene_weight.copy_(-w)
    logits_neg_on_flip = heads.forward_tissue(1.0 - mbs, present)
    assert torch.allclose(logits_raw, logits_neg_on_flip, atol=1e-5)


def test_gather_flat_region_filters_unobserved_edges() -> None:
    from mbs.training.cascade_assign import build_cascade_assignment
    from mbs.training.cascade_loop import make_synthetic_cascade_tables
    from mbs.training.flat_region_features import build_flat_region_gene_index

    tables = make_synthetic_cascade_tables(seed=7)
    assignment = build_cascade_assignment(
        locus_index=tables["locus_index"],
        locus_region_edges=tables["locus_region_edges"],
        regions=tables["regions"],
        genes=tables["genes"],
        gene_allocation="explicit_only",
    )
    index = build_flat_region_gene_index(assignment, allow_other_gene=True)
    beta = np.asarray(tables["betas"][0], dtype=np.float32)
    beta[0] = np.nan
    feats, genes = gather_flat_region_features(beta_row=beta, index=index)
    assert feats.shape[0] == genes.shape[0]
    assert feats.shape[0] < index.n_edges
    assert np.all(feats[:, -1] > 0.5)


def test_flat_region_feature_mode_ablation_zeros_role_block() -> None:
    dim = flat_region_input_dim()
    feats = np.zeros((2, dim), dtype=np.float32)
    feats[:, 1] = 1.0
    out = apply_flat_region_feature_mode(feats, "m_only")
    assert float(out[:, 1].sum()) == 0.0
