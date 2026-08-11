"""Unit tests for study-grouped auto partition (Milestone 5d)."""

from __future__ import annotations

from mbs.evaluation.splits import assert_no_study_leakage, partition_studies_by_sample_count


def test_partition_studies_by_sample_count_no_leakage() -> None:
    samples: list[dict[str, str]] = []
    for study, n in [("GSE_A", 50), ("GSE_B", 40), ("GSE_C", 30), ("GSE_D", 20), ("GSE_E", 10)]:
        samples.extend(
            {"sample_id": f"{study}_{i}", "study_id": study, "platform": "HM450"} for i in range(n)
        )
    split = partition_studies_by_sample_count(samples, seed=7, split_id="auto-test-v1")
    assert_no_study_leakage(split)
    roles = {
        "train": set(split["train_studies"]),
        "validation": set(split["validation_studies"]),
        "external_test": set(split["external_test_studies"]),
    }
    assert roles["train"] and roles["validation"] and roles["external_test"]
    assert len(roles["train"] | roles["validation"] | roles["external_test"]) == 5
    # Deterministic for fixed seed.
    split2 = partition_studies_by_sample_count(samples, seed=7, split_id="auto-test-v1")
    assert split["train_studies"] == split2["train_studies"]
