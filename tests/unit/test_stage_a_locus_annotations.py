"""Unit tests for Stage A static CpG annotation channel wiring.

Covers:
- cpg_context propagation into FlatRegionGeneIndex via build_flat_region_gene_index
- Channel vocabulary match between annotation constants and flat_region_features
- Duplicate (locus_id, gene_id) assertion fires on injected duplicate
- New feature modes (m_context, obs_only, anno_only) produce expected zero/nonzero patterns
- reg_permute_seed shuffles regulatory block while preserving M and observed
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from mbs.training.cascade_assign import build_cascade_assignment
from mbs.training.cascade_loop import make_synthetic_cascade_tables
from mbs.training.flat_region_features import (
    CPG_CONTEXTS,
    GENE_ROLES,
    PRESENCE_FLAGS,
    REGULATORY_CHANNELS,
    FlatRegionGeneIndex,
    apply_flat_region_feature_mode,
    assert_flat_region_index,
    build_flat_region_gene_index,
    flat_region_input_dim,
    gather_flat_region_features,
)


# ── helpers ───────────────────────────────────────────────────────────────────


def _make_index_with_context(ctx_label: str = "island") -> FlatRegionGeneIndex:
    """Build a FlatRegionGeneIndex with a specific cpg_context for all edges."""
    tables = make_synthetic_cascade_tables(seed=7)
    assignment = build_cascade_assignment(
        locus_index=tables["locus_index"],
        locus_region_edges=tables["locus_region_edges"],
        regions=tables["regions"],
        genes=tables["genes"],
        gene_allocation="explicit_only",
    )
    # Build a cpg_context_by_locus dict keyed on locus_id (str)
    locus_index = tables["locus_index"]
    locus_ids = locus_index["locus_id"].astype(str).tolist()
    cpg_ctx = {lid: ctx_label for lid in locus_ids}
    return build_flat_region_gene_index(
        assignment,
        locus_index=locus_index,
        cpg_context_by_locus=cpg_ctx,
        allow_other_gene=True,
    )


def _make_minimal_feats(n_edges: int) -> np.ndarray:
    """Synthetic feature matrix with known layout."""
    dim = flat_region_input_dim()
    feats = np.zeros((n_edges, dim), dtype=np.float32)
    feats[:, 0] = 0.3          # M-value
    feats[:, 1] = 1.0          # role one-hot slot 0
    feats[:, 1 + len(GENE_ROLES)] = 1.0  # ctx one-hot slot 0
    reg_start = 1 + len(GENE_ROLES) + len(CPG_CONTEXTS)
    feats[:, reg_start] = 0.5  # first regulatory channel (non-zero for permute test)
    feats[:, -1] = 1.0         # observed
    return feats


# ── vocabulary match ──────────────────────────────────────────────────────────


def test_vocabulary_sizes() -> None:
    """Channel vocabulary must match documented expected sizes."""
    assert len(GENE_ROLES) == 6, f"GENE_ROLES expected 6, got {len(GENE_ROLES)}"
    assert len(CPG_CONTEXTS) == 7, f"CPG_CONTEXTS expected 7, got {len(CPG_CONTEXTS)}"
    assert len(REGULATORY_CHANNELS) == 6, f"REGULATORY_CHANNELS expected 6, got {len(REGULATORY_CHANNELS)}"
    assert len(PRESENCE_FLAGS) == 3, f"PRESENCE_FLAGS expected 3, got {len(PRESENCE_FLAGS)}"


def test_feature_dim_formula() -> None:
    """Feature dim = 1 + n_roles + n_ctx + n_reg + n_flags + 1 (observed)."""
    expected = 1 + len(GENE_ROLES) + len(CPG_CONTEXTS) + len(REGULATORY_CHANNELS) + len(PRESENCE_FLAGS) + 1
    assert flat_region_input_dim() == expected == 24


def test_vocabulary_names_stable() -> None:
    """Vocabulary tuples must contain expected landmark values."""
    assert "promoter_core" in GENE_ROLES
    assert "gene_body" in GENE_ROLES
    assert "other_gene" in GENE_ROLES
    assert "island" in CPG_CONTEXTS
    assert "open_sea" in CPG_CONTEXTS
    assert "unknown" in CPG_CONTEXTS
    assert "promoter_like_ccre" in REGULATORY_CHANNELS
    assert "gene_role_present" in PRESENCE_FLAGS


# ── cpg_context wiring ────────────────────────────────────────────────────────


def test_cpg_context_populated_in_index() -> None:
    """When cpg_context_by_locus is provided, edge_context_id must not be all 'unknown'."""
    index = _make_index_with_context("island")
    island_id = CPG_CONTEXTS.index("island")
    # All edges should have island context
    assert np.all(index.edge_context_id == island_id), (
        "Expected all edges to be 'island' context; got: "
        + str(np.unique(index.edge_context_id))
    )
    assert np.all(index.edge_context_present), "edge_context_present should be True for island"


def test_cpg_context_unknown_when_absent() -> None:
    """Without cpg_context_by_locus, all contexts default to 'unknown'."""
    tables = make_synthetic_cascade_tables(seed=8)
    assignment = build_cascade_assignment(
        locus_index=tables["locus_index"],
        locus_region_edges=tables["locus_region_edges"],
        regions=tables["regions"],
        genes=tables["genes"],
        gene_allocation="explicit_only",
    )
    index = build_flat_region_gene_index(assignment, allow_other_gene=True)
    unknown_id = CPG_CONTEXTS.index("unknown")
    assert np.all(index.edge_context_id == unknown_id), "Default context should be 'unknown'"
    assert not np.any(index.edge_context_present), "edge_context_present should be False for unknown"


def test_graph_audit_cpg_context_counts() -> None:
    """assert_flat_region_index must report non-zero island count when context wired."""
    index = _make_index_with_context("island")
    audit = assert_flat_region_index(index)
    island_count = audit["cpg_context_counts"]["island"]
    unknown_count = audit["cpg_context_counts"]["unknown"]
    assert island_count > 0, f"Expected island count > 0, got {island_count}"
    assert unknown_count == 0, f"Expected unknown count == 0, got {unknown_count}"


# ── duplicate edge assertion ───────────────────────────────────────────────────


def test_duplicate_edge_assertion_fires() -> None:
    """assert_flat_region_index must raise on duplicate (locus_col, gene_id) pairs."""
    index = _make_index_with_context("island")
    # Inject a duplicate by repeating first row
    bad_cols = np.concatenate([index.edge_col_index, index.edge_col_index[:1]])
    bad_genes = np.concatenate([index.edge_gene_index, index.edge_gene_index[:1]])
    bad_index = FlatRegionGeneIndex(
        gene_ids=index.gene_ids,
        edge_col_index=bad_cols,
        edge_gene_index=bad_genes,
        edge_role_id=np.concatenate([index.edge_role_id, index.edge_role_id[:1]]),
        edge_context_id=np.concatenate([index.edge_context_id, index.edge_context_id[:1]]),
        edge_role_present=np.concatenate([index.edge_role_present, index.edge_role_present[:1]]),
        edge_context_present=np.concatenate([index.edge_context_present, index.edge_context_present[:1]]),
        edge_regulatory_present=np.concatenate([index.edge_regulatory_present, index.edge_regulatory_present[:1]]),
        edge_regulatory_multi_hot=np.concatenate(
            [index.edge_regulatory_multi_hot, index.edge_regulatory_multi_hot[:1]], axis=0
        ),
        n_study_loci=index.n_study_loci,
        n_other_gene_edges=index.n_other_gene_edges,
    )
    with pytest.raises(AssertionError, match="duplicate"):
        assert_flat_region_index(bad_index)


# ── new feature modes ──────────────────────────────────────────────────────────

_ROLE_START = 1
_CTX_START = _ROLE_START + len(GENE_ROLES)
_REG_START = _CTX_START + len(CPG_CONTEXTS)
_FLAGS_START = _REG_START + len(REGULATORY_CHANNELS)
_DIM = flat_region_input_dim()


def test_mode_m_context_keeps_m_and_context_zeros_role() -> None:
    n = 4
    feats = _make_minimal_feats(n)
    out = apply_flat_region_feature_mode(feats, "m_context")
    # M preserved
    np.testing.assert_array_equal(out[:, 0], feats[:, 0])
    # role one-hot zeroed
    assert out[:, _ROLE_START:_CTX_START].sum() == 0, "role block should be zero"
    # context one-hot preserved
    assert out[:, _CTX_START:_REG_START].sum() > 0, "context block should be nonzero"
    # gene_role_present flag zeroed
    assert out[:, _FLAGS_START].sum() == 0, "gene_role_present flag should be zero"
    # observed preserved
    np.testing.assert_array_equal(out[:, -1], feats[:, -1])


def test_mode_obs_only_zeros_m_and_annotation() -> None:
    n = 4
    feats = _make_minimal_feats(n)
    out = apply_flat_region_feature_mode(feats, "obs_only")
    # M zeroed
    assert out[:, 0].sum() == 0, "M should be zero in obs_only mode"
    # annotations zeroed
    assert out[:, 1:_DIM - 1].sum() == 0, "all annotation cols should be zero"
    # observed preserved
    np.testing.assert_array_equal(out[:, -1], feats[:, -1])


def test_mode_anno_only_zeros_m_preserves_role_ctx() -> None:
    n = 4
    feats = _make_minimal_feats(n)
    out = apply_flat_region_feature_mode(feats, "anno_only")
    # M zeroed
    assert out[:, 0].sum() == 0, "M should be zero in anno_only mode"
    # role block preserved
    assert out[:, _ROLE_START:_CTX_START].sum() > 0, "role block should be nonzero"
    # context block preserved
    assert out[:, _CTX_START:_REG_START].sum() > 0, "context block should be nonzero"
    # regulatory zeroed
    assert out[:, _REG_START:_FLAGS_START].sum() == 0, "regulatory block should be zero"


# ── reg_permute_seed ──────────────────────────────────────────────────────────


def test_reg_permute_shuffles_regulatory_block() -> None:
    """With reg_permute_seed, the regulatory block is shuffled; M and observed are unchanged."""
    tables = make_synthetic_cascade_tables(seed=9)
    assignment = build_cascade_assignment(
        locus_index=tables["locus_index"],
        locus_region_edges=tables["locus_region_edges"],
        regions=tables["regions"],
        genes=tables["genes"],
        gene_allocation="explicit_only",
    )
    index = build_flat_region_gene_index(assignment, allow_other_gene=True)
    # Use synthetic regulatory multi-hot to detect permutation
    n = index.n_edges
    rng_ref = np.random.default_rng(42)
    # Inject non-zero regulatory data to make permutation detectable
    reg_data = rng_ref.random((n, len(REGULATORY_CHANNELS))).astype(np.float32)
    index_with_reg = FlatRegionGeneIndex(
        gene_ids=index.gene_ids,
        edge_col_index=index.edge_col_index,
        edge_gene_index=index.edge_gene_index,
        edge_role_id=index.edge_role_id,
        edge_context_id=index.edge_context_id,
        edge_role_present=index.edge_role_present,
        edge_context_present=index.edge_context_present,
        edge_regulatory_present=index.edge_regulatory_present,
        edge_regulatory_multi_hot=reg_data,
        n_study_loci=index.n_study_loci,
        n_other_gene_edges=index.n_other_gene_edges,
    )
    beta = np.ones(index.n_study_loci, dtype=np.float32) * 0.5
    feats_base, _ = gather_flat_region_features(beta_row=beta, index=index_with_reg)
    feats_perm, _ = gather_flat_region_features(
        beta_row=beta, index=index_with_reg, reg_permute_seed=0
    )
    reg_start = _REG_START
    flags_start = _FLAGS_START
    # M and observed should be identical
    np.testing.assert_array_equal(feats_base[:, 0], feats_perm[:, 0])
    np.testing.assert_array_equal(feats_base[:, -1], feats_perm[:, -1])
    # Regulatory block should differ (with high probability for n>1 edges)
    if feats_base.shape[0] > 1:
        reg_same = np.allclose(feats_base[:, reg_start:flags_start], feats_perm[:, reg_start:flags_start])
        assert not reg_same, "Regulatory block should be shuffled with reg_permute_seed"
