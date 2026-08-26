"""Transparent gene/region-mean and elastic-net phenotype baselines (Milestone 7E)."""

from __future__ import annotations

from typing import Any, Literal

import numpy as np
from sklearn.linear_model import ElasticNet, LogisticRegression, Ridge

from mbs.evaluation.metrics import metrics_by_group, multiclass_metrics, regression_metrics

MeanKind = Literal["gene", "region"]


def presence_aware_means(
    values: np.ndarray,
    observed: np.ndarray,
    group_index: np.ndarray,
    *,
    n_groups: int,
    empty_fill: float = 0.5,
) -> np.ndarray:
    """Mean of observed values per group; empty groups → ``empty_fill``.

    ``values`` / ``observed`` are ``[n_samples, n_loci]``; ``group_index`` is
    ``[n_loci]`` mapping each locus column to a group id in ``[0, n_groups)``.
    """
    vals = np.asarray(values, dtype=np.float64)
    obs = np.asarray(observed, dtype=bool)
    groups = np.asarray(group_index, dtype=np.int64)
    if vals.shape != obs.shape:
        raise ValueError("values and observed shape mismatch")
    if groups.shape[0] != vals.shape[1]:
        raise ValueError("group_index length must match n_loci")
    n_samples = vals.shape[0]
    out = np.full((n_samples, int(n_groups)), float(empty_fill), dtype=np.float64)
    for g in range(int(n_groups)):
        cols = np.where(groups == g)[0]
        if cols.size == 0:
            continue
        block = vals[:, cols]
        mask = obs[:, cols]
        counts = mask.sum(axis=1)
        sums = np.where(mask, block, 0.0).sum(axis=1)
        has = counts > 0
        out[has, g] = sums[has] / counts[has]
    return out.astype(np.float32, copy=False)


def fit_linear_multitask(
    x_train: np.ndarray,
    *,
    age: np.ndarray | None,
    age_mask: np.ndarray | None,
    tissue: np.ndarray | None,
    tissue_mask: np.ndarray | None,
    sex: np.ndarray | None,
    sex_mask: np.ndarray | None,
) -> dict[str, Any]:
    """Fit Ridge (age) / LogReg (tissue, sex) on fixed features."""
    models: dict[str, Any] = {}
    x = np.asarray(x_train, dtype=np.float64)
    if age is not None and age_mask is not None and bool(np.asarray(age_mask).any()):
        m = np.asarray(age_mask, dtype=bool)
        models["age"] = Ridge(alpha=1.0).fit(x[m], np.asarray(age, dtype=np.float64)[m])
    if tissue is not None and tissue_mask is not None and bool(np.asarray(tissue_mask).any()):
        m = np.asarray(tissue_mask, dtype=bool)
        models["tissue"] = LogisticRegression(max_iter=500).fit(
            x[m], np.asarray(tissue, dtype=np.int64)[m]
        )
    if sex is not None and sex_mask is not None and bool(np.asarray(sex_mask).any()):
        m = np.asarray(sex_mask, dtype=bool)
        models["sex"] = LogisticRegression(max_iter=500).fit(
            x[m], np.asarray(sex, dtype=np.int64)[m]
        )
    return models


def predict_linear_multitask(
    models: dict[str, Any],
    x: np.ndarray,
) -> dict[str, np.ndarray]:
    """Predict age / tissue / sex when the corresponding model was fit."""
    x_arr = np.asarray(x, dtype=np.float64)
    out: dict[str, np.ndarray] = {}
    if "age" in models:
        out["age"] = np.asarray(models["age"].predict(x_arr), dtype=np.float64)
    if "tissue" in models:
        out["tissue"] = np.asarray(models["tissue"].predict(x_arr), dtype=np.int64)
    if "sex" in models:
        out["sex"] = np.asarray(models["sex"].predict(x_arr), dtype=np.int64)
    return out


def evaluate_multitask_predictions(
    *,
    preds: dict[str, np.ndarray],
    age: np.ndarray | None,
    age_mask: np.ndarray | None,
    tissue: np.ndarray | None,
    tissue_mask: np.ndarray | None,
    sex: np.ndarray | None,
    sex_mask: np.ndarray | None,
    study_ids: np.ndarray | None = None,
    platforms: np.ndarray | None = None,
) -> dict[str, Any]:
    """Holdout metrics for transparent / late-fusion linear heads."""
    metrics: dict[str, Any] = {}
    if "age" in preds and age is not None and age_mask is not None:
        m = np.asarray(age_mask, dtype=bool)
        if m.any():
            age_m = regression_metrics(np.asarray(age)[m], preds["age"][m])
            metrics["age"] = age_m
            if study_ids is not None:
                metrics["age_by_study"] = metrics_by_group(
                    np.asarray(age)[m],
                    preds["age"][m],
                    np.asarray(study_ids)[m],
                    task="regression",
                )
    if "tissue" in preds and tissue is not None and tissue_mask is not None:
        m = np.asarray(tissue_mask, dtype=bool)
        if m.any():
            yt = np.asarray(tissue, dtype=np.int64)[m]
            yp = preds["tissue"][m]
            metrics["tissue"] = {
                k: v
                for k, v in multiclass_metrics(yt, yp).items()
                if k in {"macro_f1", "balanced_accuracy", "accuracy"}
            }
            if study_ids is not None:
                metrics["tissue_by_study"] = metrics_by_group(
                    yt, yp, np.asarray(study_ids)[m], task="multiclass"
                )
            if platforms is not None:
                metrics["tissue_by_platform"] = metrics_by_group(
                    yt, yp, np.asarray(platforms)[m], task="multiclass"
                )
    if "sex" in preds and sex is not None and sex_mask is not None:
        m = np.asarray(sex_mask, dtype=bool)
        if m.any():
            yt = np.asarray(sex, dtype=np.int64)[m]
            yp = preds["sex"][m]
            metrics["sex"] = {
                k: v
                for k, v in multiclass_metrics(yt, yp).items()
                if k in {"macro_f1", "balanced_accuracy", "accuracy"}
            }
    return metrics


def fit_elasticnet_phenotype(
    x_train: np.ndarray,
    y_train: np.ndarray,
    *,
    task: Literal["regression", "multiclass"] = "regression",
    alpha: float = 0.1,
    l1_ratio: float = 0.5,
) -> Any:
    """Elastic-net on fixed features (transparent T-enet arm)."""
    x = np.asarray(x_train, dtype=np.float64)
    y = np.asarray(y_train)
    if task == "regression":
        model = ElasticNet(alpha=alpha, l1_ratio=l1_ratio, max_iter=2000)
        model.fit(x, y.astype(np.float64))
        return model
    # Multiclass via one-vs-rest logistic with elastic-net (saga).
    model = LogisticRegression(
        penalty="elasticnet",
        solver="saga",
        l1_ratio=l1_ratio,
        C=1.0 / max(alpha, 1e-6),
        max_iter=2000,
    )
    model.fit(x, y.astype(np.int64))
    return model


def run_mean_baseline(
    *,
    x_train: np.ndarray,
    x_test: np.ndarray,
    age_train: np.ndarray | None,
    age_mask_train: np.ndarray | None,
    tissue_train: np.ndarray | None,
    tissue_mask_train: np.ndarray | None,
    sex_train: np.ndarray | None,
    sex_mask_train: np.ndarray | None,
    age_test: np.ndarray | None,
    age_mask_test: np.ndarray | None,
    tissue_test: np.ndarray | None,
    tissue_mask_test: np.ndarray | None,
    sex_test: np.ndarray | None,
    sex_mask_test: np.ndarray | None,
    study_ids_test: np.ndarray | None = None,
    platforms_test: np.ndarray | None = None,
    kind: MeanKind = "gene",
) -> dict[str, Any]:
    """Fit + evaluate a presence-aware mean feature baseline."""
    models = fit_linear_multitask(
        x_train,
        age=age_train,
        age_mask=age_mask_train,
        tissue=tissue_train,
        tissue_mask=tissue_mask_train,
        sex=sex_train,
        sex_mask=sex_mask_train,
    )
    preds = predict_linear_multitask(models, x_test)
    metrics = evaluate_multitask_predictions(
        preds=preds,
        age=age_test,
        age_mask=age_mask_test,
        tissue=tissue_test,
        tissue_mask=tissue_mask_test,
        sex=sex_test,
        sex_mask=sex_mask_test,
        study_ids=study_ids_test,
        platforms=platforms_test,
    )
    return {"kind": kind, "metrics": metrics, "n_features": int(x_train.shape[1])}


def run_elasticnet_baseline(
    *,
    x_train: np.ndarray,
    x_test: np.ndarray,
    age_train: np.ndarray | None,
    age_mask_train: np.ndarray | None,
    tissue_train: np.ndarray | None,
    tissue_mask_train: np.ndarray | None,
    age_test: np.ndarray | None,
    age_mask_test: np.ndarray | None,
    tissue_test: np.ndarray | None,
    tissue_mask_test: np.ndarray | None,
    study_ids_test: np.ndarray | None = None,
    platforms_test: np.ndarray | None = None,
    alpha: float = 0.1,
    l1_ratio: float = 0.5,
) -> dict[str, Any]:
    """Elastic-net on gene/region mean features for age + tissue."""
    metrics: dict[str, Any] = {}
    if (
        age_train is not None
        and age_mask_train is not None
        and bool(np.asarray(age_mask_train).any())
    ):
        m = np.asarray(age_mask_train, dtype=bool)
        model = fit_elasticnet_phenotype(
            x_train[m],
            np.asarray(age_train)[m],
            task="regression",
            alpha=alpha,
            l1_ratio=l1_ratio,
        )
        if age_test is not None and age_mask_test is not None:
            mt = np.asarray(age_mask_test, dtype=bool)
            if mt.any():
                pred = model.predict(x_test[mt])
                metrics["age"] = regression_metrics(np.asarray(age_test)[mt], pred)
                if study_ids_test is not None:
                    metrics["age_by_study"] = metrics_by_group(
                        np.asarray(age_test)[mt],
                        pred,
                        np.asarray(study_ids_test)[mt],
                        task="regression",
                    )
    if (
        tissue_train is not None
        and tissue_mask_train is not None
        and bool(np.asarray(tissue_mask_train).any())
    ):
        m = np.asarray(tissue_mask_train, dtype=bool)
        model = fit_elasticnet_phenotype(
            x_train[m],
            np.asarray(tissue_train)[m],
            task="multiclass",
            alpha=alpha,
            l1_ratio=l1_ratio,
        )
        if tissue_test is not None and tissue_mask_test is not None:
            mt = np.asarray(tissue_mask_test, dtype=bool)
            if mt.any():
                pred = model.predict(x_test[mt])
                yt = np.asarray(tissue_test, dtype=np.int64)[mt]
                metrics["tissue"] = {
                    k: v
                    for k, v in multiclass_metrics(yt, pred).items()
                    if k in {"macro_f1", "balanced_accuracy", "accuracy"}
                }
                if study_ids_test is not None:
                    metrics["tissue_by_study"] = metrics_by_group(
                        yt, pred, np.asarray(study_ids_test)[mt], task="multiclass"
                    )
                if platforms_test is not None:
                    metrics["tissue_by_platform"] = metrics_by_group(
                        yt, pred, np.asarray(platforms_test)[mt], task="multiclass"
                    )
    return {"kind": "elasticnet", "metrics": metrics, "n_features": int(x_train.shape[1])}
