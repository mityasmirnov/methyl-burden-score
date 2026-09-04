"""Unit tests for CPU typed-RBS pooling (Milestone 7G), tiny synthetic arrays."""

from __future__ import annotations

import numpy as np
import pandas as pd

from mbs.training.typed_rbs_pooling import (
    apply_standardizer,
    build_layout,
    evaluate_arm,
    features_passthrough,
    features_typed,
    features_untyped,
    fit_standardizer,
    pool_gene_role,
    shuffle_region_types,
    typed_pool_promotion_gate,
)


def _tiny_region_df() -> pd.DataFrame:
    # Gene A: two promoter_core + one gene_body; Gene B: one three_prime.
    return pd.DataFrame(
        {
            "region_id": ["A:pc1", "A:pc2", "A:gb", "B:tp"],
            "gene_id": ["A", "A", "A", "B"],
            "region_type": ["promoter_core", "promoter_core", "gene_body", "three_prime"],
            "column_index": [0, 1, 2, 3],
        }
    )


def test_pool_gene_role_presence_aware() -> None:
    df = _tiny_region_df()
    layout = build_layout(df)
    assert layout.genes == ("A", "B")
    rbs = np.array([[0.9, 0.7, 0.2, 0.5], [0.6, 0.7, 0.2, 0.5]], dtype=np.float32)
    present = np.array([[1, 1, 1, 1], [1, 0, 0, 0]], dtype=np.uint8)
    pooled = pool_gene_role(rbs, present, layout)
    pc = layout.role_names.index("promoter_core")
    gb = layout.role_names.index("gene_body")
    tp = layout.role_names.index("three_prime")

    # Gene A promoter_core, sample 0: cols 0,1 present → max .4, mean .3, count 2.
    assert np.isclose(pooled["max"][0, 0, pc], 0.4)
    assert np.isclose(pooled["mean"][0, 0, pc], 0.3)
    assert pooled["count"][0, 0, pc] == 2
    assert pooled["present"][0, 0, pc] == 1.0
    # Sample 1: only col 0 present → max/mean = 0.1, count 1.
    assert np.isclose(pooled["max"][1, 0, pc], 0.1)
    assert pooled["count"][1, 0, pc] == 1
    # Gene A gene_body absent for sample 1 → neutral 0, present flag 0.
    assert pooled["present"][1, 0, gb] == 0.0
    assert pooled["max"][1, 0, gb] == 0.0
    # Gene B three_prime present for sample 0 (centered 0.0), absent sample 1.
    assert pooled["present"][0, 1, tp] == 1.0
    assert pooled["present"][1, 1, tp] == 0.0


def test_feature_shapes() -> None:
    df = _tiny_region_df()
    rbs = np.random.default_rng(0).random((5, 4)).astype(np.float32)
    present = np.ones((5, 4), dtype=np.uint8)
    x0, n0 = features_untyped(rbs, present, df)
    # R0: 1 role × 2 genes × 4 stats = 8 features.
    assert x0.shape == (5, 8)
    assert len(n0) == 8
    x1, n1 = features_typed(rbs, present, df, stats=("max",))
    # R1: (max + present + count) × 5 roles × 2 genes = 30 features.
    assert x1.shape == (5, 30)
    assert len(n1) == 30
    x3, _ = features_typed(rbs, present, df, stats=("max", "mean"))
    assert x3.shape == (5, 40)  # (max, mean, present, count) × 5 × 2


def test_passthrough_presence_aware_construction() -> None:
    df = _tiny_region_df()
    rbs = np.array([[0.9, 0.7, 0.2, 0.5]], dtype=np.float32)
    present = np.array([[1, 0, 1, 0]], dtype=np.uint8)
    feats, names = features_passthrough(rbs, present, df)
    # x = (rbs - 0.5) * present ; absent → 0.
    assert np.allclose(feats[0], [0.4, 0.0, -0.3, 0.0])
    assert len(names) == 4


def test_standardizer_train_only() -> None:
    x = np.array([[0.0, 5.0], [2.0, 5.0], [4.0, 5.0]], dtype=np.float32)
    mean, std = fit_standardizer(x)
    assert np.isclose(mean[0], 2.0)
    # Zero-variance column gets std 1 (no divide-by-zero).
    assert std[1] == 1.0
    z = apply_standardizer(x, (mean, std))
    assert np.isclose(z[:, 0].mean(), 0.0)


def test_shuffle_region_types_stats() -> None:
    df = _tiny_region_df()
    shuffled, stats = shuffle_region_types(df, seed=0)
    assert set(shuffled["region_type"]) == set(df["region_type"])  # same multiset overall
    assert 0.0 <= stats["frac_genes_altered"] <= 1.0
    assert 0.0 <= stats["frac_columns_changed"] <= 1.0
    # Single-region gene B can never change; gene A has a repeated role so the
    # fraction of altered columns is bounded below 1.
    assert stats["frac_columns_changed"] < 1.0


def test_promotion_gate_promote() -> None:
    gate = typed_pool_promotion_gate(
        age_mae_r0=10.0,
        age_mae_typed=8.0,
        age_r2_r0=0.5,
        age_r2_typed=0.5,
        tissue_f1_r0=0.50,
        tissue_f1_typed=0.50,
        sex_auroc_r0=None,
        sex_auroc_typed=None,
        shuffle_age_mae=9.5,  # ≥1 y worse than typed
    )
    assert gate["promote"] is True


def test_promotion_gate_blockers() -> None:
    base = {
        "age_mae_r0": 10.0,
        "age_mae_typed": 8.0,
        "age_r2_r0": 0.5,
        "age_r2_typed": 0.5,
        "tissue_f1_r0": 0.5,
        "tissue_f1_typed": 0.5,
        "sex_auroc_r0": 0.9,
        "sex_auroc_typed": 0.9,
        "shuffle_age_mae": 9.5,
    }
    # No age gain.
    g = typed_pool_promotion_gate(**{**base, "age_mae_typed": 9.7})
    assert g["promote"] is False and g["age_improves"] is False
    # Tissue drops too much.
    g = typed_pool_promotion_gate(**{**base, "tissue_f1_typed": 0.4})
    assert g["promote"] is False and g["tissue_ok"] is False
    # Sex AUROC drops too much.
    g = typed_pool_promotion_gate(**{**base, "sex_auroc_typed": 0.8})
    assert g["promote"] is False and g["sex_ok"] is False
    # Shuffle is not clearly worse than typed (≥1 y required).
    g = typed_pool_promotion_gate(**{**base, "shuffle_age_mae": 8.5})
    assert g["promote"] is False and g["shuffle_worse_than_typed"] is False
    # R2 path alone can promote even without MAE gain (shuffle still clearly worse).
    g = typed_pool_promotion_gate(
        **{**base, "age_mae_typed": 9.7, "age_r2_typed": 0.56, "shuffle_age_mae": 11.0}
    )
    assert g["promote"] is True and g["age_improves"] is True


def test_evaluate_arm_end_to_end_age_only() -> None:
    rng = np.random.default_rng(1)
    df = _tiny_region_df()
    n = 24
    rbs = rng.random((n, 4)).astype(np.float32)
    present = (rng.random((n, 4)) > 0.2).astype(np.uint8)
    feats, _ = features_typed(rbs, present, df, stats=("max", "mean"))
    age = (rng.random(n) * 40 + 20).astype(np.float64)
    train_idx = np.arange(0, 16)
    test_idx = np.arange(16, n)
    falsemask = np.zeros(n, dtype=bool)
    arrays = {
        "age": age,
        "age_mask": np.ones(n, dtype=bool),
        "tissue": np.zeros(n, dtype=np.int64),
        "tissue_mask": falsemask,
        "sex": np.zeros(n, dtype=np.int64),
        "sex_mask": falsemask,
        "study_ids": np.array([f"s{i % 3}" for i in range(n)], dtype=object),
    }
    result = evaluate_arm(
        name="R3",
        x_train=feats[train_idx],
        x_test=feats[test_idx],
        pheno_train={k: v[train_idx] for k, v in arrays.items()},
        pheno_test={k: v[test_idx] for k, v in arrays.items()},
    )
    assert result["arm"] == "R3"
    assert "mae" in result["metrics"]["age"]
