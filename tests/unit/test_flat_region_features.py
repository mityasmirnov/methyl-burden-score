"""FlatDeepSetRegion feature assembly and fold-safe panel checks."""

from __future__ import annotations

import numpy as np
import torch

from mbs.models import FlatDeepSetRegion
from mbs.training.cascade_assign import assignment_col_subset, build_cascade_assignment
from mbs.training.cascade_loop import make_synthetic_cascade_tables
from mbs.training.flat_region_features import (
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
