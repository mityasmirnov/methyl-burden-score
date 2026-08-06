"""Unit tests for evaluation metrics and study-grouped splits."""

from __future__ import annotations

import numpy as np
import pytest

from mbs.evaluation import (
    assert_no_study_leakage,
    binary_auroc_auprc,
    build_study_grouped_split,
    multiclass_metrics,
    regression_metrics,
)


def test_regression_metrics_perfect() -> None:
    y = np.array([1.0, 2.0, 3.0])
    out = regression_metrics(y, y)
    assert out["mae"] == 0.0
    assert out["rmse"] == 0.0


def test_binary_auroc_separable() -> None:
    y = np.array([0, 0, 1, 1])
    s = np.array([0.1, 0.2, 0.8, 0.9])
    out = binary_auroc_auprc(y, s)
    assert out["auroc"] == pytest.approx(1.0)
    assert out["auprc"] == pytest.approx(1.0)


def test_multiclass_metrics_balanced() -> None:
    y = np.array([0, 1, 2, 0, 1, 2])
    p = np.array([0, 1, 2, 0, 1, 2])
    out = multiclass_metrics(y, p, n_classes=3)
    assert out["macro_f1"] == pytest.approx(1.0)
    assert out["balanced_accuracy"] == pytest.approx(1.0)
    assert out["confusion_matrix"] == [[2, 0, 0], [0, 2, 0], [0, 0, 2]]


def test_study_grouped_split_rejects_leakage() -> None:
    samples = [
        {"sample_id": "a1", "study_id": "G1"},
        {"sample_id": "a2", "study_id": "G1"},
        {"sample_id": "b1", "study_id": "G2"},
    ]
    with pytest.raises(ValueError, match="leakage"):
        build_study_grouped_split(
            samples,
            train_studies=["G1"],
            validation_studies=["G1"],
            external_test_studies=["G2"],
        )


def test_study_grouped_split_ok() -> None:
    samples = [
        {"sample_id": "a1", "study_id": "G1", "platform": "HM450"},
        {"sample_id": "b1", "study_id": "G2", "platform": "HM450"},
        {"sample_id": "c1", "study_id": "G3", "platform": "EPIC"},
    ]
    split = build_study_grouped_split(
        samples,
        train_studies=["G1"],
        validation_studies=["G2"],
        external_test_studies=["G3"],
    )
    assert_no_study_leakage(split)
    assert split["train_sample_ids"] == ["a1"]
    assert split["external_test_sample_ids"] == ["c1"]
    assert split["mode"] == "study_grouped"
