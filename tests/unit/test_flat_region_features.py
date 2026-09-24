"""FlatDeepSetRegion feature assembly and fold-safe panel checks."""

from __future__ import annotations

import numpy as np
import torch

from mbs.models import FlatDeepSetRegion
from mbs.training.cascade_assign import assignment_col_subset, build_cascade_assignment
from mbs.training.cascade_loop import make_synthetic_cascade_tables
from mbs.training.flat_region_features import (
    build_flat_region_base_features,
    build_flat_region_gene_index,
    flat_region_input_dim,
    gather_flat_region_features,
)
from mbs.training.fold_safe_panel import expand_panel_columns, stability_select_columns


def test_flat_region_input_dim_and_gather() -> None:
    tables = make_synthetic_cascade_tables(seed=1)
    assignment = build_cascade_assignment(
        locus_index=tables["locus_index"],
        locus_region_edges=tables["locus_region_edges"],
        regions=tables["regions"],
        genes=tables["genes"],
        gene_allocation="explicit_only",
    )
    index = build_flat_region_gene_index(assignment, allow_other_gene=True)
    dim = flat_region_input_dim()
    beta_row = np.asarray(tables["betas"][0], dtype=np.float32)
    feats, cpg_to_gene = gather_flat_region_features(beta_row=beta_row, index=index)
    assert feats.shape[1] == dim
    assert cpg_to_gene.shape[0] == feats.shape[0]
    model = FlatDeepSetRegion(dim, phi_hidden_dim=8, rho_hidden_dim=4, phi_layers=1, rho_layers=1)
    x = torch.from_numpy(feats)
    g = torch.from_numpy(cpg_to_gene)
    out = model(x, g, index.n_genes)
    assert out["mbs"].shape == (index.n_genes,)


def test_flat_region_static_block_appends_and_survives_ablation() -> None:
    """CpGPT-style static embedding: appended after the base 24 cols, untouched
    by feature_mode ablation (m_only zeroes annotation but not static)."""
    tables = make_synthetic_cascade_tables(seed=3)
    assignment = build_cascade_assignment(
        locus_index=tables["locus_index"],
        locus_region_edges=tables["locus_region_edges"],
        regions=tables["regions"],
        genes=tables["genes"],
        gene_allocation="explicit_only",
    )
    index = build_flat_region_gene_index(assignment, allow_other_gene=True)
    base_dim = flat_region_input_dim()
    static_dim = 5
    rng = np.random.default_rng(0)
    n_cols = int(tables["betas"].shape[1])
    static_by_col = rng.normal(size=(n_cols, static_dim)).astype(np.float32)
    edge_static_block = static_by_col[index.edge_col_index]

    base = build_flat_region_base_features(
        index, feature_mode="m_only", static_block=edge_static_block
    )
    assert base.shape == (index.n_edges, base_dim + static_dim)
    # Annotation block (cols 1..base_dim-2) zeroed by m_only; static survives.
    assert np.allclose(base[:, 1 : base_dim - 1], 0.0)
    assert np.allclose(base[:, base_dim:], edge_static_block)

    beta_row = np.asarray(tables["betas"][0], dtype=np.float32)
    feats, cpg_to_gene = gather_flat_region_features(
        beta_row=beta_row,
        index=index,
        base_features=base,
        feature_mode="m_only",
        static_dim=static_dim,
    )
    assert feats.shape[1] == base_dim + static_dim
    assert cpg_to_gene.shape[0] == feats.shape[0]
    # Observed flag still lands at base_dim - 1, not swallowed by the static block.
    assert set(np.unique(feats[:, base_dim - 1]).tolist()) <= {0.0, 1.0}
    # Static columns for observed edges match the source block exactly.
    obs = np.isfinite(beta_row[index.edge_col_index])
    assert np.allclose(feats[:, base_dim:], edge_static_block[obs])


def test_flat_region_static_dim_without_base_features_raises() -> None:
    tables = make_synthetic_cascade_tables(seed=4)
    assignment = build_cascade_assignment(
        locus_index=tables["locus_index"],
        locus_region_edges=tables["locus_region_edges"],
        regions=tables["regions"],
        genes=tables["genes"],
        gene_allocation="explicit_only",
    )
    index = build_flat_region_gene_index(assignment, allow_other_gene=True)
    beta_row = np.asarray(tables["betas"][0], dtype=np.float32)
    try:
        gather_flat_region_features(beta_row=beta_row, index=index, static_dim=3)
    except ValueError as error:
        assert "static_dim" in str(error)
    else:
        raise AssertionError("expected ValueError for static_dim without base_features")


def test_assignment_col_subset_and_panel_expand() -> None:
    tables = make_synthetic_cascade_tables(seed=2)
    assignment = build_cascade_assignment(
        locus_index=tables["locus_index"],
        locus_region_edges=tables["locus_region_edges"],
        regions=tables["regions"],
        genes=tables["genes"],
    )
    cols = np.asarray([0], dtype=np.int64)
    sub = assignment_col_subset(assignment, cols)
    assert sub.edge_col_index.size >= 0
    expanded = expand_panel_columns(cols, assignment)
    assert expanded.size >= 1


def test_stability_select_columns_smoke() -> None:
    rng = np.random.default_rng(0)
    x = rng.normal(size=(40, 20)).astype(np.float32)
    y = rng.integers(0, 3, size=40)
    picked, meta = stability_select_columns(x, y, max_seeds=5, n_inner_folds=2)
    assert picked.size >= 1
    assert meta["n_runs"] >= 0
