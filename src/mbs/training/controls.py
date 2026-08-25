"""Negative-control transforms (static-only, coverage-only, label permutation)."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Sequence
from random import Random
from typing import Any

import numpy as np
from numpy.typing import NDArray
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.preprocessing import OneHotEncoder

from mbs.evaluation.metrics import (
    binary_auroc_auprc,
    masked_multilabel_auroc_auprc,
    multiclass_metrics,
    regression_metrics,
)
from mbs.training.phenotypes import MultilabelMaps, SamplePhenotype


def apply_feature_control(
    features: NDArray[np.float32],
    *,
    mode: str,
    include_m_value: bool = True,
    include_robust_z: bool = False,
) -> NDArray[np.float32]:
    """Zero methylation and/or static channels.

    Layout: ``beta, [M], [z], static..., static_present, [norm_present]``.
    Coverage-only keeps trailing present flags (``static_present`` and optional
    ``norm_present``).
    """
    out = np.asarray(features, dtype=np.float32).copy()
    if out.ndim != 2 or out.shape[0] == 0:
        return out
    n_methyl = 1
    if include_m_value:
        n_methyl += 1
    if include_robust_z:
        n_methyl += 1
    n_flags = 1 + (1 if include_robust_z else 0)
    if mode in {"none", "off", ""}:
        return out
    if mode == "static_only":
        out[:, :n_methyl] = 0.0
        return out
    if mode == "coverage_only":
        # Keep only trailing present flag column(s).
        out[:, : out.shape[1] - n_flags] = 0.0
        return out
    raise ValueError(f"unknown feature control: {mode}")


def permute_labels_within_study(
    phenotypes: Sequence[SamplePhenotype],
    *,
    seed: int,
) -> list[SamplePhenotype]:
    """Shuffle class_index / age / sex within study strata."""
    rng = Random(seed)  # noqa: S311
    by_study: dict[str, list[int]] = defaultdict(list)
    for i, ph in enumerate(phenotypes):
        by_study[str(ph.study_id or ph.sample_id)].append(i)
    ages = [p.age for p in phenotypes]
    classes = [p.class_index for p in phenotypes]
    sexes = [p.sex_class_index for p in phenotypes]
    for idxs in by_study.values():
        order = list(range(len(idxs)))
        rng.shuffle(order)
        src_ages = [ages[i] for i in idxs]
        src_cls = [classes[i] for i in idxs]
        src_sex = [sexes[i] for i in idxs]
        for k, i in enumerate(idxs):
            ages[i] = src_ages[order[k]]
            classes[i] = src_cls[order[k]]
            sexes[i] = src_sex[order[k]]
    return [
        SamplePhenotype(
            sample_id=ph.sample_id,
            cell_type=ph.cell_type,
            donor_id=ph.donor_id,
            title=ph.title,
            class_index=classes[i],
            study_id=ph.study_id,
            age=ages[i],
            platform=ph.platform,
            age_mask=ph.age_mask,
            tissue_mask=ph.tissue_mask,
            sex_mask=ph.sex_mask,
            sex_class_index=sexes[i],
        )
        for i, ph in enumerate(phenotypes)
    ]


def _meta_design(
    study_ids: Sequence[str],
    platforms: Sequence[str | None],
    tissues: Sequence[str | None],
    *,
    encoder: OneHotEncoder | None = None,
) -> tuple[Any, OneHotEncoder]:
    cats = np.column_stack(
        [
            np.asarray(study_ids, dtype=object),
            np.asarray([p or "unknown" for p in platforms], dtype=object),
            np.asarray([t or "unknown" for t in tissues], dtype=object),
        ]
    )
    if encoder is None:
        enc = OneHotEncoder(handle_unknown="ignore")
        x = enc.fit_transform(cats)
        return x, enc
    return encoder.transform(cats), encoder


def fit_metadata_only(
    *,
    study_ids: Sequence[str],
    platforms: Sequence[str | None],
    tissues: Sequence[str | None],
    y: np.ndarray,
    task: str,
    eval_study_ids: Sequence[str] | None = None,
    eval_platforms: Sequence[str | None] | None = None,
    eval_tissues: Sequence[str | None] | None = None,
    eval_y: np.ndarray | None = None,
) -> dict[str, float]:
    """Linear/logistic ceiling from study/platform/tissue one-hots (not the encoder).

    When eval_* are provided, fit on the train fold and score the holdout
    (OOF confounding ceiling). Otherwise fit+score on the same rows (legacy).
    """
    x_train, enc = _meta_design(study_ids, platforms, tissues)
    y_arr = np.asarray(y)
    if task == "regression":
        model = Ridge(alpha=1.0).fit(x_train, y_arr)
        if eval_y is None:
            pred = model.predict(x_train)
            return regression_metrics(y_arr, pred)
        x_eval, _ = _meta_design(
            eval_study_ids or (),
            eval_platforms or (),
            eval_tissues or (),
            encoder=enc,
        )
        return regression_metrics(np.asarray(eval_y), model.predict(x_eval))
    model = LogisticRegression(max_iter=200).fit(x_train, y_arr)
    if eval_y is None:
        pred = model.predict(x_train)
        return {
            k: v
            for k, v in multiclass_metrics(y_arr, pred).items()
            if k in {"macro_f1", "balanced_accuracy"}
        }
    x_eval, _ = _meta_design(
        eval_study_ids or (),
        eval_platforms or (),
        eval_tissues or (),
        encoder=enc,
    )
    pred = model.predict(x_eval)
    return {
        k: v
        for k, v in multiclass_metrics(np.asarray(eval_y), pred).items()
        if k in {"macro_f1", "balanced_accuracy"}
    }


def _ph_meta(ph: SamplePhenotype) -> tuple[str, str | None, str | None]:
    return (
        str(ph.study_id or ph.sample_id),
        ph.platform,
        ph.cell_type if ph.tissue_mask else None,
    )


def evaluate_metadata_only_ceiling(
    *,
    train: Sequence[SamplePhenotype],
    eval_sets: dict[str, Sequence[SamplePhenotype]],
    disease_maps: MultilabelMaps | None = None,
    cancer_maps: MultilabelMaps | None = None,
) -> dict[str, Any]:
    """Train-fit metadata ceiling; score each named eval fold (val/test).

    This is a sidecar confounding report — it does not replace neural training.
    """
    out: dict[str, Any] = {"protocol": "fit_train_score_holdout"}
    # Age
    train_age = [p for p in train if p.age_mask and p.age is not None]
    if len(train_age) >= 2:
        for fold_name, fold_ph in eval_sets.items():
            eval_age = [p for p in fold_ph if p.age_mask and p.age is not None]
            if len(eval_age) < 1:
                continue
            train_ages = [float(a) for p in train_age if (a := p.age) is not None]
            eval_ages = [float(a) for p in eval_age if (a := p.age) is not None]
            metrics = fit_metadata_only(
                study_ids=[_ph_meta(p)[0] for p in train_age],
                platforms=[_ph_meta(p)[1] for p in train_age],
                tissues=[_ph_meta(p)[2] for p in train_age],
                y=np.array(train_ages),
                task="regression",
                eval_study_ids=[_ph_meta(p)[0] for p in eval_age],
                eval_platforms=[_ph_meta(p)[1] for p in eval_age],
                eval_tissues=[_ph_meta(p)[2] for p in eval_age],
                eval_y=np.array(eval_ages),
            )
            out.setdefault(fold_name, {})["age"] = metrics

    # Tissue
    train_tis = [p for p in train if p.tissue_mask]
    if len(train_tis) >= 2 and len({p.class_index for p in train_tis}) >= 2:
        for fold_name, fold_ph in eval_sets.items():
            eval_tis = [p for p in fold_ph if p.tissue_mask]
            if len(eval_tis) < 1:
                continue
            metrics = fit_metadata_only(
                study_ids=[_ph_meta(p)[0] for p in train_tis],
                platforms=[_ph_meta(p)[1] for p in train_tis],
                tissues=[_ph_meta(p)[2] for p in train_tis],
                y=np.array([p.class_index for p in train_tis]),
                task="multiclass",
                eval_study_ids=[_ph_meta(p)[0] for p in eval_tis],
                eval_platforms=[_ph_meta(p)[1] for p in eval_tis],
                eval_tissues=[_ph_meta(p)[2] for p in eval_tis],
                eval_y=np.array([p.class_index for p in eval_tis]),
            )
            out.setdefault(fold_name, {})["tissue"] = metrics

    # Sex
    train_sex = [p for p in train if p.sex_mask]
    if len(train_sex) >= 2 and len({p.sex_class_index for p in train_sex}) >= 2:
        for fold_name, fold_ph in eval_sets.items():
            eval_sex = [p for p in fold_ph if p.sex_mask]
            if len(eval_sex) < 1:
                continue
            metrics = fit_metadata_only(
                study_ids=[_ph_meta(p)[0] for p in train_sex],
                platforms=[_ph_meta(p)[1] for p in train_sex],
                tissues=[_ph_meta(p)[2] for p in train_sex],
                y=np.array([p.sex_class_index for p in train_sex]),
                task="multiclass",
                eval_study_ids=[_ph_meta(p)[0] for p in eval_sex],
                eval_platforms=[_ph_meta(p)[1] for p in eval_sex],
                eval_tissues=[_ph_meta(p)[2] for p in eval_sex],
                eval_y=np.array([p.sex_class_index for p in eval_sex]),
            )
            out.setdefault(fold_name, {})["sex"] = metrics

    # Disease / cancer multilabel: one-vs-rest logistic per label with both classes.
    for task_name, maps in (("disease", disease_maps), ("cancer", cancer_maps)):
        if maps is None or not maps.label_names:
            continue
        for fold_name, fold_ph in eval_sets.items():
            label_metrics: dict[str, Any] = {}
            for li, lab in enumerate(maps.label_names):
                train_ids = [
                    p.sample_id
                    for p in train
                    if maps.masks.get(p.sample_id) is not None and bool(maps.masks[p.sample_id][li])
                ]
                eval_ids = [
                    p.sample_id
                    for p in fold_ph
                    if maps.masks.get(p.sample_id) is not None and bool(maps.masks[p.sample_id][li])
                ]
                if len(train_ids) < 2 or len(eval_ids) < 1:
                    continue
                y_tr = np.array([float(maps.targets[s][li]) for s in train_ids])
                y_ev = np.array([float(maps.targets[s][li]) for s in eval_ids])
                if len(np.unique(y_tr)) < 2 or len(np.unique(y_ev)) < 2:
                    continue
                ph_by = {p.sample_id: p for p in [*train, *fold_ph]}
                metrics = fit_metadata_only(
                    study_ids=[_ph_meta(ph_by[s])[0] for s in train_ids],
                    platforms=[_ph_meta(ph_by[s])[1] for s in train_ids],
                    tissues=[_ph_meta(ph_by[s])[2] for s in train_ids],
                    y=y_tr.astype(np.int64),
                    task="multiclass",
                    eval_study_ids=[_ph_meta(ph_by[s])[0] for s in eval_ids],
                    eval_platforms=[_ph_meta(ph_by[s])[1] for s in eval_ids],
                    eval_tissues=[_ph_meta(ph_by[s])[2] for s in eval_ids],
                    eval_y=y_ev.astype(np.int64),
                )
                # Also try probability AUROC when both classes exist.
                try:
                    x_tr, enc = _meta_design(
                        [_ph_meta(ph_by[s])[0] for s in train_ids],
                        [_ph_meta(ph_by[s])[1] for s in train_ids],
                        [_ph_meta(ph_by[s])[2] for s in train_ids],
                    )
                    clf = LogisticRegression(max_iter=200).fit(x_tr, y_tr.astype(np.int64))
                    x_ev, _ = _meta_design(
                        [_ph_meta(ph_by[s])[0] for s in eval_ids],
                        [_ph_meta(ph_by[s])[1] for s in eval_ids],
                        [_ph_meta(ph_by[s])[2] for s in eval_ids],
                        encoder=enc,
                    )
                    scores = clf.predict_proba(x_ev)[:, 1]
                    metrics.update(binary_auroc_auprc(y_ev.astype(np.int64), scores))
                except (ValueError, IndexError):
                    pass
                label_metrics[str(lab)] = metrics
            if label_metrics:
                # Macro over labels that produced metrics.
                aurocs = [float(v["auroc"]) for v in label_metrics.values() if "auroc" in v]
                fold_summary: dict[str, Any] = {"n_labels_scored": len(label_metrics)}
                if aurocs:
                    fold_summary["macro_auroc"] = float(np.mean(aurocs))
                out.setdefault(fold_name, {})[task_name] = fold_summary
                # Keep multilabel helper available for future dense tables.
                _ = masked_multilabel_auroc_auprc
    return out
