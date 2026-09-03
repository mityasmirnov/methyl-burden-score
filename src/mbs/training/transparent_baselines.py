"""Transparent gene/region-mean and elastic-net phenotype baselines (Milestone 7E)."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any, Literal

import numpy as np
from sklearn.decomposition import PCA
from sklearn.linear_model import ElasticNet, LogisticRegression, Ridge, SGDClassifier

from mbs.evaluation.metrics import metrics_by_group, multiclass_metrics, regression_metrics

MeanKind = Literal["gene", "region"]
TissueSolver = Literal["logistic", "balanced_logistic", "sgd_ovr"]


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


def _fit_tissue_classifier(
    x: np.ndarray,
    y: np.ndarray,
    *,
    tissue_solver: TissueSolver = "logistic",
) -> Any:
    """Multiclass tissue head for transparent / late-fusion baselines."""
    if tissue_solver == "logistic":
        return LogisticRegression(max_iter=1000).fit(x, y)
    if tissue_solver == "balanced_logistic":
        return LogisticRegression(max_iter=1000, class_weight="balanced").fit(x, y)
    if tissue_solver == "sgd_ovr":
        return SGDClassifier(
            loss="log_loss",
            class_weight="balanced",
            max_iter=1000,
            random_state=0,
        ).fit(x, y)
    raise ValueError(f"unknown tissue_solver: {tissue_solver}")


def fit_linear_multitask(
    x_train: np.ndarray,
    *,
    age: np.ndarray | None,
    age_mask: np.ndarray | None,
    tissue: np.ndarray | None,
    tissue_mask: np.ndarray | None,
    sex: np.ndarray | None,
    sex_mask: np.ndarray | None,
    tissue_solver: TissueSolver = "logistic",
    fusion_pca_components: int | None = None,
) -> dict[str, Any]:
    """Fit Ridge (age) / LogReg (tissue, sex) on fixed features."""
    models: dict[str, Any] = {}
    x = np.asarray(x_train, dtype=np.float64)
    x_tissue = x
    if fusion_pca_components is not None:
        n_comp = min(int(fusion_pca_components), x.shape[0], x.shape[1])
        if n_comp < 1:
            raise ValueError("fusion_pca_components must be positive")
        pca = PCA(n_components=n_comp, random_state=0)
        pca.fit(x)
        x_tissue = pca.transform(x)
        models["_pca"] = pca
    if age is not None and age_mask is not None and bool(np.asarray(age_mask).any()):
        m = np.asarray(age_mask, dtype=bool)
        models["age"] = Ridge(alpha=1.0).fit(x[m], np.asarray(age, dtype=np.float64)[m])
    if tissue is not None and tissue_mask is not None and bool(np.asarray(tissue_mask).any()):
        m = np.asarray(tissue_mask, dtype=bool)
        models["tissue"] = _fit_tissue_classifier(
            x_tissue[m],
            np.asarray(tissue, dtype=np.int64)[m],
            tissue_solver=tissue_solver,
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
    x_tissue = x_arr
    if "_pca" in models:
        x_tissue = models["_pca"].transform(x_arr)
    out: dict[str, np.ndarray] = {}
    if "age" in models:
        out["age"] = np.asarray(models["age"].predict(x_arr), dtype=np.float64)
    if "tissue" in models:
        out["tissue"] = np.asarray(models["tissue"].predict(x_tissue), dtype=np.int64)
        if hasattr(models["tissue"], "predict_proba"):
            out["tissue_proba"] = np.asarray(
                models["tissue"].predict_proba(x_tissue), dtype=np.float64
            )
            out["tissue_classes"] = np.asarray(models["tissue"].classes_, dtype=np.int64)
    if "sex" in models:
        out["sex"] = np.asarray(models["sex"].predict(x_arr), dtype=np.int64)
        if hasattr(models["sex"], "predict_proba"):
            out["sex_proba"] = np.asarray(models["sex"].predict_proba(x_arr), dtype=np.float64)
            out["sex_classes"] = np.asarray(models["sex"].classes_, dtype=np.int64)
    return out


def _positive_class_scores(
    proba: np.ndarray,
    classes: np.ndarray,
    *,
    positive: int = 1,
) -> np.ndarray | None:
    """Column of ``predict_proba`` for ``positive`` class, if present."""
    classes_a = np.asarray(classes, dtype=np.int64)
    proba_a = np.asarray(proba, dtype=np.float64)
    matches = np.where(classes_a == int(positive))[0]
    if matches.size == 0:
        if proba_a.ndim == 2 and proba_a.shape[1] == 2:
            return proba_a[:, 1]
        return None
    return proba_a[:, int(matches[0])]


def tissue_one_vs_rest_auroc(
    y_true: np.ndarray,
    proba: np.ndarray,
    classes: np.ndarray,
    *,
    top_n: int = 5,
    class_names: list[str] | None = None,
) -> list[dict[str, Any]]:
    """Top-N tissue one-vs-rest AUROC curves from multiclass probabilities."""
    from sklearn.metrics import roc_auc_score, roc_curve  # noqa: PLC0415

    yt = np.asarray(y_true, dtype=np.int64)
    proba_a = np.asarray(proba, dtype=np.float64)
    classes_a = np.asarray(classes, dtype=np.int64)
    class_to_col = {int(c): i for i, c in enumerate(classes_a.tolist())}
    counts: dict[int, int] = {}
    for yi in yt.tolist():
        counts[int(yi)] = counts.get(int(yi), 0) + 1
    ranked = sorted(counts.items(), key=lambda kv: -kv[1])[: int(top_n)]
    curves: list[dict[str, Any]] = []
    for cls, _n in ranked:
        if cls not in class_to_col:
            continue
        binary = (yt == cls).astype(np.int64)
        if binary.min() == binary.max():
            continue
        scores = proba_a[:, class_to_col[cls]]
        fpr, tpr, _ = roc_curve(binary, scores)
        label = (
            class_names[cls]
            if class_names is not None and 0 <= cls < len(class_names)
            else str(cls)
        )
        # Downsample for report JSON size.
        if fpr.size > 80:
            idx = np.linspace(0, fpr.size - 1, 80).astype(np.int64)
            fpr_s, tpr_s = fpr[idx], tpr[idx]
        else:
            fpr_s, tpr_s = fpr, tpr
        curves.append(
            {
                "label": label,
                "class_index": int(cls),
                "auroc": float(roc_auc_score(binary, scores)),
                "fpr": fpr_s.tolist(),
                "tpr": tpr_s.tolist(),
            }
        )
    return curves


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
    tissue_class_names: list[str] | None = None,
    tissue_valid_classes: Iterable[int] | None = None,
) -> dict[str, Any]:
    """Holdout metrics for transparent / late-fusion linear heads.

    ``tissue_valid_classes`` should be the set of tissue class indices seen
    during training for this fold; classes never seen in train cannot be
    predicted correctly by construction and are excluded from macro-F1 /
    balanced-accuracy (see ``multiclass_metrics``).
    """
    from mbs.evaluation.metrics import binary_auroc_auprc  # noqa: PLC0415

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
                for k, v in multiclass_metrics(yt, yp, valid_classes=tissue_valid_classes).items()
                if k
                in {
                    "macro_f1",
                    "balanced_accuracy",
                    "accuracy",
                    "n_classes_scored",
                    "excluded_zero_shot_test_counts",
                }
            }
            if (
                "tissue_proba" in preds
                and "tissue_classes" in preds
                and preds["tissue_proba"].shape[0] == len(preds["tissue"])
            ):
                proba_m = preds["tissue_proba"][m]
                curves = tissue_one_vs_rest_auroc(
                    yt,
                    proba_m,
                    preds["tissue_classes"],
                    class_names=tissue_class_names,
                )
                if curves:
                    metrics["tissue_roc"] = curves
                    # Mean AUROC over only the top-N most frequent classes
                    # (see tissue_one_vs_rest_auroc); NOT the same class set as
                    # macro_f1/balanced_accuracy above, which cover every class
                    # present in y_true. Do not compare this value directly
                    # against macro_f1 as if they shared a denominator.
                    metrics["tissue"]["top5_ovr_auroc"] = float(
                        np.mean([c["auroc"] for c in curves])
                    )
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
            sex_metrics: dict[str, Any] = {
                k: v
                for k, v in multiclass_metrics(yt, yp).items()
                if k in {"macro_f1", "balanced_accuracy", "accuracy"}
            }
            if (
                "sex_proba" in preds
                and "sex_classes" in preds
                and preds["sex_proba"].shape[0] == len(preds["sex"])
            ):
                scores = _positive_class_scores(
                    preds["sex_proba"][m], preds["sex_classes"], positive=1
                )
                if scores is not None and len(np.unique(yt)) >= 2:
                    try:
                        from sklearn.metrics import roc_curve  # noqa: PLC0415

                        roc = binary_auroc_auprc(yt, scores)
                        sex_metrics.update(roc)
                        fpr, tpr, _ = roc_curve(yt, scores)
                        if fpr.size > 80:
                            idx = np.linspace(0, fpr.size - 1, 80).astype(np.int64)
                            fpr, tpr = fpr[idx], tpr[idx]
                        sex_metrics["fpr"] = fpr.tolist()
                        sex_metrics["tpr"] = tpr.tolist()
                    except ValueError:
                        pass
            metrics["sex"] = sex_metrics
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
    tissue_class_names: list[str] | None = None,
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
    tissue_valid_classes = None
    if tissue_train is not None and tissue_mask_train is not None:
        tm = np.asarray(tissue_mask_train, dtype=bool)
        if tm.any():
            tissue_valid_classes = set(np.asarray(tissue_train, dtype=np.int64)[tm].tolist())
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
        tissue_class_names=tissue_class_names,
        tissue_valid_classes=tissue_valid_classes,
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
                tissue_valid_classes = set(np.asarray(tissue_train, dtype=np.int64)[m].tolist())
                metrics["tissue"] = {
                    k: v
                    for k, v in multiclass_metrics(
                        yt, pred, valid_classes=tissue_valid_classes
                    ).items()
                    if k
                    in {
                        "macro_f1",
                        "balanced_accuracy",
                        "accuracy",
                        "n_classes_scored",
                        "excluded_zero_shot_test_counts",
                    }
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


def run_elasticnet_multitask(
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
    tissue_class_names: list[str] | None = None,
    alpha: float = 0.1,
    l1_ratio: float = 0.5,
) -> dict[str, Any]:
    """Fold-fit elastic-net heads on frozen features (age + tissue + sex)."""
    models: dict[str, Any] = {}
    if (
        age_train is not None
        and age_mask_train is not None
        and bool(np.asarray(age_mask_train).any())
    ):
        m = np.asarray(age_mask_train, dtype=bool)
        models["age"] = fit_elasticnet_phenotype(
            x_train[m],
            np.asarray(age_train)[m],
            task="regression",
            alpha=alpha,
            l1_ratio=l1_ratio,
        )
    if (
        tissue_train is not None
        and tissue_mask_train is not None
        and bool(np.asarray(tissue_mask_train).any())
    ):
        m = np.asarray(tissue_mask_train, dtype=bool)
        models["tissue"] = fit_elasticnet_phenotype(
            x_train[m],
            np.asarray(tissue_train)[m],
            task="multiclass",
            alpha=alpha,
            l1_ratio=l1_ratio,
        )
    if (
        sex_train is not None
        and sex_mask_train is not None
        and bool(np.asarray(sex_mask_train).any())
    ):
        m = np.asarray(sex_mask_train, dtype=bool)
        models["sex"] = fit_elasticnet_phenotype(
            x_train[m],
            np.asarray(sex_train)[m],
            task="multiclass",
            alpha=alpha,
            l1_ratio=l1_ratio,
        )
    preds = predict_linear_multitask(models, x_test)
    tissue_valid_classes = None
    if tissue_train is not None and tissue_mask_train is not None:
        tm = np.asarray(tissue_mask_train, dtype=bool)
        if tm.any():
            tissue_valid_classes = set(np.asarray(tissue_train, dtype=np.int64)[tm].tolist())
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
        tissue_class_names=tissue_class_names,
        tissue_valid_classes=tissue_valid_classes,
    )
    return {
        "kind": "elasticnet",
        "metrics": metrics,
        "n_features": int(np.asarray(x_train).shape[1]),
        "alpha": float(alpha),
        "l1_ratio": float(l1_ratio),
    }
