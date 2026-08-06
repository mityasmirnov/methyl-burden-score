"""Study-grouped train / validation / external_test split manifests."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from random import Random
from typing import Any


@dataclass(frozen=True, slots=True)
class SampleSplitRow:
    sample_id: str
    study_id: str
    role: str  # train | validation | external_test
    platform: str | None = None
    group_id: str | None = None


def _roles_by_study(
    study_ids: Iterable[str],
    *,
    train_studies: Sequence[str],
    validation_studies: Sequence[str],
    external_test_studies: Sequence[str],
) -> dict[str, str]:
    train = {str(x) for x in train_studies}
    val = {str(x) for x in validation_studies}
    test = {str(x) for x in external_test_studies}
    overlap_tv = train & val
    overlap_te = train & test
    overlap_ve = val & test
    if overlap_tv or overlap_te or overlap_ve:
        raise ValueError(
            "study leakage across roles: "
            f"train∩val={sorted(overlap_tv)} train∩test={sorted(overlap_te)} "
            f"val∩test={sorted(overlap_ve)}"
        )
    mapping: dict[str, str] = {}
    for sid in train:
        mapping[sid] = "train"
    for sid in val:
        mapping[sid] = "validation"
    for sid in test:
        mapping[sid] = "external_test"
    unknown = {str(s) for s in study_ids} - set(mapping)
    if unknown:
        raise ValueError(f"studies missing from split roles: {sorted(unknown)}")
    return mapping


def build_study_grouped_split(
    samples: Sequence[dict[str, Any]],
    *,
    train_studies: Sequence[str],
    validation_studies: Sequence[str],
    external_test_studies: Sequence[str] | None = None,
    split_id: str = "study-grouped-v1",
) -> dict[str, Any]:
    """Build a leakage-safe study-grouped split manifest.

    ``samples`` entries must include ``sample_id`` and ``study_id``; optional
    ``platform`` and ``group_id`` (defaults to study_id).
    """
    external = list(external_test_studies or [])
    study_ids = [str(s["study_id"]) for s in samples]
    role_by_study = _roles_by_study(
        study_ids,
        train_studies=train_studies,
        validation_studies=validation_studies,
        external_test_studies=external,
    )
    rows: list[dict[str, Any]] = []
    for sample in samples:
        study_id = str(sample["study_id"])
        sample_id = str(sample["sample_id"])
        role = role_by_study[study_id]
        group_id = str(sample.get("group_id") or study_id)
        rows.append(
            {
                "sample_id": sample_id,
                "study_id": study_id,
                "platform": sample.get("platform"),
                "group_id": group_id,
                "role": role,
            }
        )
    # Enforce: no sample role conflicts; no study spanning roles (already ensured).
    by_role: dict[str, list[str]] = {"train": [], "validation": [], "external_test": []}
    for row in rows:
        by_role[row["role"]].append(row["sample_id"])
    return {
        "split_id": split_id,
        "mode": "study_grouped",
        "train_studies": sorted({str(x) for x in train_studies}),
        "validation_studies": sorted({str(x) for x in validation_studies}),
        "external_test_studies": sorted({str(x) for x in external}),
        "samples": rows,
        "train_sample_ids": by_role["train"],
        "validation_sample_ids": by_role["validation"],
        "external_test_sample_ids": by_role["external_test"],
    }


def assert_no_study_leakage(split: dict[str, Any]) -> None:
    """Raise if any study_id appears under more than one role."""
    role_by_study: dict[str, str] = {}
    for row in split.get("samples", []):
        study_id = str(row["study_id"])
        role = str(row["role"])
        prev = role_by_study.get(study_id)
        if prev is not None and prev != role:
            raise ValueError(f"study {study_id} appears in both {prev} and {role}")
        role_by_study[study_id] = role


def partition_studies_by_sample_count(
    samples: Sequence[dict[str, Any]],
    *,
    seed: int = 42,
    train_fraction: float = 0.7,
    val_fraction: float = 0.15,
    split_id: str = "study-grouped-auto-v1",
) -> dict[str, Any]:
    """Seeded study-grouped split targeting sample-count fractions.

    Studies are sorted by descending sample count (then study_id), shuffled with
    a seeded RNG after a stable sort key, then greedily assigned to train → val
    → test until target sample budgets are met. Every study appears in exactly
    one role.
    """
    if train_fraction <= 0 or val_fraction < 0 or train_fraction + val_fraction >= 1.0:
        raise ValueError("need train_fraction > 0, val_fraction >= 0, sum < 1")
    counts: dict[str, int] = {}
    for sample in samples:
        study_id = str(sample["study_id"])
        counts[study_id] = counts.get(study_id, 0) + 1
    if not counts:
        raise ValueError("no samples to partition")
    total = sum(counts.values())
    train_budget = int(total * train_fraction)
    val_budget = int(total * val_fraction)
    # Ensure non-empty roles when enough studies exist.
    studies = sorted(counts.keys(), key=lambda s: (-counts[s], s))
    rng = Random(seed)
    rng.shuffle(studies)

    train: list[str] = []
    val: list[str] = []
    test: list[str] = []
    train_n = val_n = 0
    for study in studies:
        n = counts[study]
        if train_n < train_budget or not train:
            train.append(study)
            train_n += n
        elif val_n < val_budget or not val:
            val.append(study)
            val_n += n
        else:
            test.append(study)
    if not test and len(train) > 1:
        moved = train.pop()
        test.append(moved)
        train_n -= counts[moved]
    if not val and len(train) > 1:
        moved = train.pop()
        val.append(moved)
        train_n -= counts[moved]
    if not train or not val or not test:
        raise ValueError(
            f"auto split produced empty role(s): "
            f"train={len(train)} val={len(val)} test={len(test)} "
            f"(n_studies={len(counts)})"
        )
    return build_study_grouped_split(
        samples,
        train_studies=train,
        validation_studies=val,
        external_test_studies=test,
        split_id=split_id,
    )
