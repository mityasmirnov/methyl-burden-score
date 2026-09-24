"""CpGPT/DNA-LM static-embedding plumbing for the cascade encoder.

Mirrors the flat_region contract: static embeddings are appended as **trailing**
columns after the M-value, so column 0 stays the M-value and the historical
``[..., 1]`` layout is preserved byte-for-byte when no static block is given.
Covers both training and scoring, which share the single ``_forward_batch``
funnel.
"""

from __future__ import annotations

import numpy as np
import pytest
import torch

from mbs.models import CascadeDeepSet
from mbs.training.cascade_assign import build_cascade_assignment
from mbs.training.cascade_loop import (
    _dense_cpg_features_batch,
    _forward_batch,
    make_synthetic_cascade_tables,
)


def test_no_static_block_is_unchanged() -> None:
    """static_edge_block=None must keep the historical [batch, n_edges, 1] layout."""
    betas = np.array([[0.1, 0.9, np.nan], [0.5, 0.2, 0.7]], dtype=np.float32)
    feats = _dense_cpg_features_batch(betas)
    assert feats.shape == (2, 3, 1)
    # explicit None and a width-0 block agree with the default
    assert np.array_equal(feats, _dense_cpg_features_batch(betas, static_edge_block=None))
    empty = np.zeros((3, 0), dtype=np.float32)
    assert np.array_equal(feats, _dense_cpg_features_batch(betas, static_edge_block=empty))


def test_static_block_appends_trailing_columns() -> None:
    """M-value stays at column 0; static dims land after it and are shared per edge."""
    betas = np.array([[0.1, 0.9, 0.4], [0.5, 0.2, 0.7]], dtype=np.float32)
    base = _dense_cpg_features_batch(betas)
    static = np.arange(3 * 4, dtype=np.float32).reshape(3, 4)  # [n_edges, static_dim]
    feats = _dense_cpg_features_batch(betas, static_edge_block=static)

    assert feats.shape == (2, 3, 5)  # 1 M-value + 4 static
    np.testing.assert_allclose(feats[..., :1], base)  # column 0 unchanged
    for b in range(betas.shape[0]):  # static is per-edge, broadcast across samples
        np.testing.assert_allclose(feats[b, :, 1:], static)


def test_static_block_shape_mismatch_fails_closed() -> None:
    betas = np.zeros((2, 3), dtype=np.float32)
    with pytest.raises(ValueError, match="static_edge_block must be"):
        _dense_cpg_features_batch(betas, static_edge_block=np.zeros((7, 4), dtype=np.float32))


def test_forward_batch_accepts_static_block_end_to_end() -> None:
    """The real forward pass runs with static dims and matches the widened input_dim."""
    tables = make_synthetic_cascade_tables(seed=5)
    assignment = build_cascade_assignment(
        locus_index=tables["locus_index"],
        locus_region_edges=tables["locus_region_edges"],
        regions=tables["regions"],
        genes=tables["genes"],
    )
    betas = np.asarray(tables["betas"], dtype=np.float32)[:2]
    n_edges = int(assignment.edge_col_index.size)
    static_dim = 6
    rng = np.random.default_rng(0)
    static = rng.normal(size=(n_edges, static_dim)).astype(np.float32)

    device = torch.device("cpu")
    n_region_types = max(len(assignment.region_types), 1)
    model = CascadeDeepSet(1 + static_dim, n_region_types, cpg_hidden_dim=8, region_hidden_dim=4)
    out = _forward_batch(
        model, assignment, betas, device=device, static_edge_block=static
    )
    assert out["mbs"].shape[0] == betas.shape[0]
    assert torch.isfinite(out["mbs"]).all()

    # And the width-1 model still works with no static block (regression guard).
    plain = CascadeDeepSet(1, n_region_types, cpg_hidden_dim=8, region_hidden_dim=4)
    out_plain = _forward_batch(plain, assignment, betas, device=device)
    assert out_plain["mbs"].shape == out["mbs"].shape
    assert torch.isfinite(out_plain["mbs"]).all()
