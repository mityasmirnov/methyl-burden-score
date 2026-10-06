"""Within-gene epoch sampler + sparse gather helpers (§3.5)."""

from __future__ import annotations

import numpy as np

from mbs.training.flat_region_features import (
    FlatRegionGeneIndex,
    gather_flat_region_features,
)
from mbs.training.sparse_betas import SparseBetasView
from mbs.training.within_gene_sampler import WithinGeneEpochSampler


def _toy_index() -> FlatRegionGeneIndex:
    from mbs.training.flat_region_features import REGULATORY_CHANNELS  # noqa: PLC0415

    n_reg = len(REGULATORY_CHANNELS)
    gene_ids = ["g0", "g1", "g2"]
    edge_gene = np.asarray([0, 0, 1, 1, 1, 1, 1, 2], dtype=np.int64)
    n = edge_gene.size
    return FlatRegionGeneIndex(
        gene_ids=gene_ids,
        edge_col_index=np.arange(n, dtype=np.int64),
        edge_gene_index=edge_gene,
        edge_role_id=np.zeros(n, dtype=np.int64),
        edge_context_id=np.zeros(n, dtype=np.int64),
        edge_role_present=np.ones(n, dtype=bool),
        edge_context_present=np.ones(n, dtype=bool),
        edge_regulatory_present=np.zeros(n, dtype=bool),
        edge_regulatory_multi_hot=np.zeros((n, n_reg), dtype=np.float32),
        n_study_loci=n,
        n_other_gene_edges=0,
    )


def test_epoch_sampler_resamples_and_tracks_exposure() -> None:
    index = _toy_index()
    sampler = WithinGeneEpochSampler(full_index=index, max_cpgs_per_gene=2, seed=0)
    a = sampler.sample_epoch(1)
    b = sampler.sample_epoch(2)
    assert a.n_genes == 3 and b.n_genes == 3
    assert a.n_edges == 2 + 2 + 1
    # Gene 1 (5 edges) should not keep the exact same 2 edges for every epoch forever.
    # Across two epochs at least one edge differs or exposure accumulates.
    exp = sampler.exposure_report()
    assert exp["n_epochs"] == 2
    assert exp["times_selected_sum"] == a.n_edges + b.n_edges
    sel = np.asarray(exp["per_edge"]["times_selected"])
    assert int(sel.sum()) == a.n_edges + b.n_edges
    assert np.all(np.asarray(exp["per_edge"]["times_available"]) == 2)


def test_gather_accepts_beta_edge() -> None:
    index = _toy_index()
    beta_edge = np.linspace(0.1, 0.9, index.n_edges).astype(np.float32)
    feats, genes = gather_flat_region_features(index=index, beta_edge=beta_edge)
    assert feats.shape[0] == index.n_edges
    assert genes.shape[0] == feats.shape[0]


class _FakeHandle:
    def __init__(self, arr: np.ndarray) -> None:
        self.arr = arr

    def __getitem__(self, key: object) -> np.ndarray:  # noqa: ANN401
        return self.arr[key]


def test_sparse_betas_view_no_universe_preload() -> None:
    mat = np.arange(20, dtype=np.float32).reshape(4, 5)
    view = SparseBetasView(handle=_FakeHandle(mat), n_cols=5)
    block = view.gather_rows_cols(np.asarray([0, 2]), np.asarray([1, 3]))
    assert block.shape == (2, 2)
    assert float(block[0, 0]) == float(mat[0, 1])
    edge = view.beta_edge_values(1, np.asarray([0, 0, 4], dtype=np.int64))
    assert edge.shape == (3,)
    assert float(edge[2]) == float(mat[1, 4])
    assert view.stats.n_calls >= 2
