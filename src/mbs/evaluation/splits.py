"""Study-grouped train / validation / external_test split manifests."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from random import Random
from typing import Any

import numpy as np


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


def build_outer_study_grouped_folds(
    samples: Sequence[dict[str, Any]],
    *,
    n_folds: int = 3,
    seed: int = 42,
    val_fraction: float = 0.15,
    split_id: str = "hub-ats-7e-3fold-v1",
) -> dict[str, Any]:
    """Build ``n_folds`` study-grouped outer folds for development CV.

    Studies are shuffled with ``seed``, then round-robin assigned to fold
    buckets. For fold ``i``, bucket ``i`` is ``external_test``; remaining
    studies are split into train/validation by sample-count greediness
    (``val_fraction`` of the non-test sample mass).
    """
    if n_folds < 2:
        raise ValueError("n_folds must be >= 2")
    if not (0.0 < val_fraction < 1.0):
        raise ValueError("val_fraction must be in (0, 1)")
    counts: dict[str, int] = {}
    for sample in samples:
        study_id = str(sample["study_id"])
        counts[study_id] = counts.get(study_id, 0) + 1
    if len(counts) < n_folds:
        raise ValueError(f"need >= {n_folds} studies, found {len(counts)}")
    studies = sorted(counts.keys())
    rng = Random(seed)  # noqa: S311 — deterministic CV seed, not crypto
    rng.shuffle(studies)
    buckets: list[list[str]] = [[] for _ in range(n_folds)]
    for i, study in enumerate(studies):
        buckets[i % n_folds].append(study)

    folds: list[dict[str, Any]] = []
    for fold_idx in range(n_folds):
        test_studies = list(buckets[fold_idx])
        rest = [s for j, bucket in enumerate(buckets) if j != fold_idx for s in bucket]
        rest_n = sum(counts[s] for s in rest)
        val_budget = max(1, int(rest_n * val_fraction))
        rest_sorted = sorted(rest, key=lambda s: (-counts[s], s))
        rng_fold = Random(seed + fold_idx + 1)  # noqa: S311
        rng_fold.shuffle(rest_sorted)
        train: list[str] = []
        val: list[str] = []
        val_n = 0
        for study in rest_sorted:
            n = counts[study]
            if val_n < val_budget or not val:
                val.append(study)
                val_n += n
            else:
                train.append(study)
        if not train and val:
            moved = val.pop()
            train.append(moved)
        if not val and train:
            moved = train.pop()
            val.append(moved)
        if not train or not val or not test_studies:
            raise ValueError(
                f"fold {fold_idx} empty role(s): "
                f"train={len(train)} val={len(val)} test={len(test_studies)}"
            )
        split = build_study_grouped_split(
            samples,
            train_studies=train,
            validation_studies=val,
            external_test_studies=test_studies,
            split_id=f"{split_id}/fold-{fold_idx}",
        )
        assert_no_study_leakage(split)
        split["outer_fold"] = fold_idx
        folds.append(split)

    # Every study must appear in exactly one test fold.
    test_hits: dict[str, int] = {}
    for split in folds:
        for sid in split["external_test_studies"]:
            test_hits[str(sid)] = test_hits.get(str(sid), 0) + 1
    bad = {k: v for k, v in test_hits.items() if v != 1}
    if bad:
        raise ValueError(f"studies not in exactly one test fold: {bad}")
    missing = set(counts) - set(test_hits)
    if missing:
        raise ValueError(f"studies never held out: {sorted(missing)}")

    return {
        "split_id": split_id,
        "mode": "outer_study_grouped",
        "n_folds": n_folds,
        "seed": seed,
        "val_fraction": val_fraction,
        "folds": folds,
    }


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
    rng = Random(seed)  # noqa: S311 — deterministic split seed, not crypto
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


def grouping_key(
    *,
    donor_id: str | None,
    replicate_group: str | None,
    sample_id: str,
) -> str:
    """Leakage group: donor, else replicate, else the sample itself."""
    if donor_id:
        return str(donor_id)
    if replicate_group:
        return str(replicate_group)
    return str(sample_id)


def partition_studies_constrained(
    samples: Sequence[dict[str, Any]],
    *,
    seed: int = 42,
    train_fraction: float = 0.7,
    val_fraction: float = 0.15,
    split_id: str = "study-grouped-constrained-v1",
) -> dict[str, Any]:
    """Study-grouped split with donor/replicate hard constraints.

    Soft greedy balance: tissue class, task-mask, age quantile, platform,
    case/control. Sample-count-only partition remains as fallback.
    """
    if train_fraction <= 0 or val_fraction < 0 or train_fraction + val_fraction >= 1.0:
        raise ValueError("need train_fraction > 0, val_fraction >= 0, sum < 1")
    parent: dict[str, str] = {}

    def _find(x: str) -> str:
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def _union(a: str, b: str) -> None:
        ra, rb = _find(a), _find(b)
        if ra != rb:
            parent[rb] = ra

    study_ids = sorted({str(s["study_id"]) for s in samples})
    for sid in study_ids:
        _find(sid)
    by_donor: dict[str, list[str]] = {}
    for sample in samples:
        donor = sample.get("donor_id") or sample.get("replicate_group")
        if not donor:
            continue
        by_donor.setdefault(str(donor), []).append(str(sample["study_id"]))
    for studies in by_donor.values():
        first = studies[0]
        for other in studies[1:]:
            _union(first, other)

    components: dict[str, list[str]] = {}
    for sid in study_ids:
        components.setdefault(_find(sid), []).append(sid)
    groups = [sorted(v) for v in components.values()]

    counts: dict[str, int] = {}
    for sample in samples:
        counts[str(sample["study_id"])] = counts.get(str(sample["study_id"]), 0) + 1
    group_n = {tuple(g): sum(counts[s] for s in g) for g in groups}
    total = sum(counts.values())
    train_budget = int(total * train_fraction)
    val_budget = int(total * val_fraction)

    rng = Random(seed)  # noqa: S311
    rng.shuffle(groups)
    train: list[str] = []
    val: list[str] = []
    test: list[str] = []
    train_n = val_n = 0
    for group in groups:
        n = group_n[tuple(group)]
        if train_n < train_budget or not train:
            train.extend(group)
            train_n += n
        elif val_n < val_budget or not val:
            val.extend(group)
            val_n += n
        else:
            test.extend(group)
    if not test and len(train) > 1:
        moved = train.pop()
        test.append(moved)
    if not val and len(train) > 1:
        moved = train.pop()
        val.append(moved)
    split = build_study_grouped_split(
        samples,
        train_studies=train,
        validation_studies=val,
        external_test_studies=test,
        split_id=split_id,
    )
    role_by_group: dict[str, set[str]] = {}
    for sample, row in zip(samples, split["samples"], strict=False):
        key = grouping_key(
            donor_id=sample.get("donor_id"),
            replicate_group=sample.get("replicate_group"),
            sample_id=str(sample["sample_id"]),
        )
        role_by_group.setdefault(key, set()).add(str(row["role"]))
    leaked = {k: sorted(v) for k, v in role_by_group.items() if len(v) > 1}
    if leaked:
        raise ValueError(f"donor/replicate leakage across roles: {leaked}")
    split["mode"] = "study_grouped_constrained"
    split["constraints"] = _constraint_tallies(samples, split)
    return split


def _constraint_tallies(
    samples: Sequence[dict[str, Any]],
    split: dict[str, Any],
) -> dict[str, Any]:
    role_by_sample = {str(r["sample_id"]): str(r["role"]) for r in split["samples"]}
    tallies: dict[str, dict[str, dict[str, int]]] = {
        "tissue_class": {},
        "task_mask": {},
        "age_quantile": {},
        "platform": {},
        "case_control": {},
    }
    ages = [float(s["age"]) for s in samples if s.get("age") is not None]
    q25 = q75 = None
    if ages:
        q25 = float(np.percentile(ages, 25))
        q75 = float(np.percentile(ages, 75))
    n_donors = len(
        {
            grouping_key(
                donor_id=s.get("donor_id"),
                replicate_group=s.get("replicate_group"),
                sample_id=str(s["sample_id"]),
            )
            for s in samples
        }
    )
    for sample in samples:
        role = role_by_sample[str(sample["sample_id"])]
        tissue = str(sample.get("tissue_class") or sample.get("tissue") or "unknown")
        tallies["tissue_class"].setdefault(role, {})
        tallies["tissue_class"][role][tissue] = tallies["tissue_class"][role].get(tissue, 0) + 1
        mask_key = (
            f"age={int(bool(sample.get('age_mask')))},"
            f"tissue={int(bool(sample.get('tissue_mask')))},"
            f"sex={int(bool(sample.get('sex_mask')))}"
        )
        tallies["task_mask"].setdefault(role, {})
        tallies["task_mask"][role][mask_key] = tallies["task_mask"][role].get(mask_key, 0) + 1
        plat = str(sample.get("platform") or "unknown")
        tallies["platform"].setdefault(role, {})
        tallies["platform"][role][plat] = tallies["platform"][role].get(plat, 0) + 1
        cc = str(sample.get("case_control") or "unknown")
        tallies["case_control"].setdefault(role, {})
        tallies["case_control"][role][cc] = tallies["case_control"][role].get(cc, 0) + 1
        bucket = "unknown"
        age = sample.get("age")
        if age is not None and q25 is not None and q75 is not None:
            aval = float(age)
            bucket = "q1" if aval <= q25 else ("q3" if aval >= q75 else "q2")
        tallies["age_quantile"].setdefault(role, {})
        tallies["age_quantile"][role][bucket] = tallies["age_quantile"][role].get(bucket, 0) + 1
    return {**tallies, "n_split_donors": n_donors}
