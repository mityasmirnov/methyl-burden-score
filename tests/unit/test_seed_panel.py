"""Unit checks for the internal_fold seed-gene panel constructor."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from mbs.training.cascade_assign import CascadeAssignment
from mbs.training.seed_panel import (
    SeedPanelArtifacts,
    build_internal_fold_seed_panel,
    gene_mask_tensor,
    matched_random_gene_panel,
    write_seed_panel,
)


def _fake_assignment() -> CascadeAssignment:
    """Tiny explicit-edge assignment: G0 has cols {0,1}, G1={2}, G2={3}, orphan=4, direct=5."""
    return CascadeAssignment(
        gene_ids=["G0", "G1", "G2"],
        region_ids=["G0:prom", "G0:body", "G1:body", "G2:body", "RBS:orphan"],
        region_type_id=np.array([0, 4, 4, 4, 5], dtype=np.int64),
        region_to_gene=np.array([0, 0, 1, 2, -1], dtype=np.int64),
        orphan_region_mask=np.array([False, False, False, False, True]),
        edge_col_index=np.array([0, 1, 2, 3, 4], dtype=np.int64),
        edge_region_index=np.array([0, 1, 2, 3, 4], dtype=np.int64),
        direct_col_index=np.array([5], dtype=np.int64),
        region_types=(
            "promoter_core",
            "promoter_proximal",
            "five_prime",
            "three_prime",
            "gene_body",
            "cgi_island",
        ),
        n_study_loci=6,
        allocated_gene_id=["G0", "G0", "G1", "G2", None],
    )


def _age_signal_data(seed: int = 0) -> dict[str, np.ndarray]:
    """12 samples, 3 studies; age is a strong linear function of col 0 only."""
    rng = np.random.default_rng(seed)
    n, n_cols = 12, 6
    x = rng.uniform(0.05, 0.95, size=(n, n_cols)).astype(np.float32)
    # col 1 is a zero-variance G0 sibling: never a stability seed, always enriched.
    x[:, 1] = 0.5
    age = (40.0 + 40.0 * x[:, 0] + rng.normal(0, 0.05, size=n)).astype(np.float64)
    study_ids = np.asarray([f"ST{i // 4}" for i in range(n)], dtype=object)
    return {"x": x, "age": age, "study_ids": study_ids}


def _build_age_only(min_genes: int, n_genes: int = 10, seed: int = 0) -> SeedPanelArtifacts:
    data = _age_signal_data(seed)
    n = data["x"].shape[0]
    zeros = np.zeros(n, dtype=np.int64)
    false_mask = np.zeros(n, dtype=bool)
    return build_internal_fold_seed_panel(
        x_train=data["x"],
        age=data["age"],
        age_mask=np.ones(n, dtype=bool),
        sex=zeros,
        sex_mask=false_mask,
        tissue=zeros,
        tissue_mask=false_mask,
        study_ids=data["study_ids"],
        assignment=_fake_assignment(),
        n_genes=n_genes,
        min_genes=min_genes,
    )


def test_gene_enrichment_includes_sibling_cpgs() -> None:
    artifacts = _build_age_only(min_genes=1)
    loci = artifacts.loci
    g0 = loci[(loci["gene_id"] == "G0") & (loci["trait"] == "age")]
    # G0 must appear (col 0 is the age signal), enriched with its sibling col 1.
    assert set(g0["locus_col"].tolist()) == {0, 1}
    seed_flags = dict(zip(g0["locus_col"], g0["is_seed_cpg"], strict=True))
    assert bool(seed_flags[0])  # col 0 is the seed
    assert not bool(seed_flags[1])  # col 1 enriched only, not selected
    # Sanity: the seed panel JSON round-trips the trait summary.
    assert artifacts.panel_json["traits"]["age"]["n_genes_actual"] >= 1


def test_undersized_panel_raises() -> None:
    # Only 3 genes exist in the assignment; require more -> must fail loudly.
    with pytest.raises(ValueError, match="min_genes"):
        _build_age_only(min_genes=5)


def test_matched_random_excludes_seed_genes() -> None:
    seed_genes = ["G0", "G1"]
    candidates = ["G0", "G1", "G2", "G3", "G4"]
    cpg_counts = {"G0": 10, "G1": 3, "G2": 9, "G3": 4, "G4": 20}
    matched = matched_random_gene_panel(
        seed_genes,
        candidate_gene_ids=candidates,
        gene_cpg_counts=cpg_counts,
        rng=np.random.default_rng(0),
    )
    assert len(matched) == len(seed_genes)
    assert set(matched).isdisjoint(seed_genes)
    assert len(set(matched)) == len(matched)  # no replacement
    # Nearest CpG count: G0(10)->G2(9), G1(3)->G3(4).
    assert matched == ["G2", "G3"]


def test_matched_random_insufficient_candidates_raises() -> None:
    with pytest.raises(ValueError, match="candidates"):
        matched_random_gene_panel(
            ["G0", "G1"],
            candidate_gene_ids=["G0", "G1", "G2"],
            gene_cpg_counts={"G2": 1},
            rng=np.random.default_rng(0),
        )


def test_gene_mask_tensor_shape_and_zeros() -> None:
    mask = gene_mask_tensor([1, 3], n_genes=5, n_outputs=2)
    assert mask.shape == (2, 5)
    assert mask.dtype == np.float32
    np.testing.assert_array_equal(mask[0], np.array([0, 1, 0, 1, 0], dtype=np.float32))
    np.testing.assert_array_equal(mask[0], mask[1])
    assert gene_mask_tensor([], n_genes=4).sum() == 0.0
    with pytest.raises(ValueError, match="out of range"):
        gene_mask_tensor([9], n_genes=4)


def test_write_seed_panel_round_trip(tmp_path: Path) -> None:
    artifacts = _build_age_only(min_genes=1)
    paths = write_seed_panel(tmp_path, artifacts)
    assert paths["panel_hash"] == artifacts.panel_hash
    assert (tmp_path / "seed_panel.json").is_file()
    assert (tmp_path / "seed_panel_gene.parquet").is_file()
    assert (tmp_path / "seed_panel_locus.parquet").is_file()
