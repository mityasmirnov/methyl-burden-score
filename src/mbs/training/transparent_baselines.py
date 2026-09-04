"""Transparent gene/region-mean and elastic-net phenotype baselines (Milestone 7E)."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any, Literal

import numpy as np
from sklearn.decomposition import PCA
from sklearn.linear_model import ElasticNet, LogisticRegression, Ridge, SGDClassifier, SGDRegressor
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score
from sklearn.preprocessing import StandardScaler

from mbs.evaluation.metrics import metrics_by_group, multiclass_metrics, regression_metrics
from mbs.training.fold_safe_panel import _study_inner_folds

MeanKind = Literal["gene", "region"]
TissueSolver = Literal["logistic", "balanced_logistic", "sgd_ovr"]

# Nested elastic-net grid (includes the fixed diagnostic α=0.1 as one candidate).
NESTED_ENET_ALPHA_GRID: tuple[float, ...] = (1e-5, 1e-4, 1e-3, 1e-2, 0.1)
NESTED_ENET_L1_GRID: tuple[float, ...] = (0.25, 0.5, 0.75)
_WIDE_FEATURE_THRESHOLD = 512


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
        return LogisticRegression(max_iter=1000, n_jobs=-1).fit(x, y)
    if tissue_solver == "balanced_logistic":
        return LogisticRegression(max_iter=1000, class_weight="balanced", n_jobs=-1).fit(x, y)
    if tissue_solver == "sgd_ovr":
        return SGDClassifier(
            loss="log_loss",
            class_weight="balanced",
            max_iter=1000,
            random_state=0,
            n_jobs=-1,
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
    max_iter: int | None = None,
) -> Any:
    """Elastic-net on fixed features (transparent T-enet arm)."""
    x = np.asarray(x_train, dtype=np.float64)
    y = np.asarray(y_train)
    n_features = int(x.shape[1]) if x.ndim == 2 else 0
    # ponytail: saga OVR on gene-MBS (~2–3k cols, many tissues) is multi-hour;
    # SGD for Stage A diagnostics and any moderately wide panel.
    wide = n_features > 512
    if task == "regression":
        iters = 400 if wide else (max_iter if max_iter is not None else 2000)
        model = ElasticNet(alpha=alpha, l1_ratio=l1_ratio, max_iter=iters)
        model.fit(x, y.astype(np.float64))
        return model
    if wide:
        # SGDClassifier.alpha ≠ LogisticRegression C=1/alpha. Scale by n_samples so
        # Stage A gene-MBS (~2–3k dims) matches short-saga tissue F1 in ~20s not hours.
        iters = max_iter if max_iter is not None else 800
        sgd_alpha = float(alpha) / max(int(x.shape[0]), 1)
        model = SGDClassifier(
            loss="log_loss",
            penalty="elasticnet",
            alpha=max(sgd_alpha, 1e-8),
            l1_ratio=float(l1_ratio),
            max_iter=iters,
            tol=1e-3,
            random_state=0,
            n_jobs=-1,
        )
        model.fit(x, y.astype(np.int64))
        return model
    # Multiclass via one-vs-rest logistic with elastic-net (saga).
    iters = max_iter if max_iter is not None else 2000
    model = LogisticRegression(
        penalty="elasticnet",
        solver="saga",
        l1_ratio=l1_ratio,
        C=1.0 / max(alpha, 1e-6),
        max_iter=iters,
        n_jobs=-1,
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


def _fit_nested_enet_model(
    x_train: np.ndarray,
    y_train: np.ndarray,
    *,
    task: Literal["regression", "multiclass"],
    alpha: float,
    l1_ratio: float,
) -> Any:
    """Fit elastic-net for nested readout; SGD on wide panels (RBS ~13k cols)."""
    x = np.asarray(x_train, dtype=np.float64)
    y = np.asarray(y_train)
    wide = int(x.shape[1]) > _WIDE_FEATURE_THRESHOLD
    if task == "regression":
        if wide:
            model = SGDRegressor(
                loss="squared_error",
                penalty="elasticnet",
                alpha=float(alpha),
                l1_ratio=float(l1_ratio),
                max_iter=200,
                tol=1e-3,
                random_state=0,
            )
            model.fit(x, y.astype(np.float64))
            return model
        model = ElasticNet(alpha=alpha, l1_ratio=l1_ratio, max_iter=2000)
        model.fit(x, y.astype(np.float64))
        return model
    # multiclass
    if wide:
        sgd_alpha = float(alpha) / max(int(x.shape[0]), 1)
        model = SGDClassifier(
            loss="log_loss",
            penalty="elasticnet",
            alpha=max(sgd_alpha, 1e-8),
            l1_ratio=float(l1_ratio),
            max_iter=200,
            tol=1e-3,
            random_state=0,
            n_jobs=-1,
        )
        model.fit(x, y.astype(np.int64))
        return model
    model = LogisticRegression(
        penalty="elasticnet",
        solver="saga",
        l1_ratio=l1_ratio,
        C=1.0 / max(alpha, 1e-6),
        max_iter=2000,
        n_jobs=-1,
    )
    model.fit(x, y.astype(np.int64))
    return model


def _nested_grids_for_width(
    n_features: int,
    *,
    alpha_grid: tuple[float, ...] = NESTED_ENET_ALPHA_GRID,
    l1_grid: tuple[float, ...] = NESTED_ENET_L1_GRID,
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Shrink hyperparam grid for very wide RBS (~13k) so nested stays practical."""
    if n_features > 2048:
        # Keep fixed diagnostic α=0.1 in the set; drop the finest alphas.
        alphas = tuple(a for a in alpha_grid if a >= 1e-3)
        return alphas or alpha_grid, l1_grid
    return alpha_grid, l1_grid


def _scale_train_test(
    x_train: np.ndarray,
    x_test: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, StandardScaler]:
    """Fit StandardScaler on train only; transform train and test."""
    scaler = StandardScaler(with_mean=True, with_std=True)
    x_tr = scaler.fit_transform(np.asarray(x_train, dtype=np.float64))
    x_te = scaler.transform(np.asarray(x_test, dtype=np.float64))
    return x_tr, x_te, scaler


def _inner_val_score(
    *,
    task: Literal["age", "tissue", "sex"],
    y_true: np.ndarray,
    y_pred: np.ndarray,
    y_proba: np.ndarray | None = None,
) -> float:
    """Higher is better inner-validation score for hyperparameter selection."""
    if task == "age":
        yt = np.asarray(y_true, dtype=np.float64)
        yp = np.asarray(y_pred, dtype=np.float64)
        mae = float(np.mean(np.abs(yt - yp)))
        return -mae
    if task == "tissue":
        return float(
            f1_score(
                np.asarray(y_true, dtype=np.int64),
                np.asarray(y_pred, dtype=np.int64),
                average="macro",
                zero_division=0,
            )
        )
    # sex
    if y_proba is not None and y_proba.ndim == 2 and y_proba.shape[1] >= 2:
        try:
            return float(roc_auc_score(np.asarray(y_true, dtype=np.int64), y_proba[:, 1]))
        except ValueError:
            pass
    return float(
        accuracy_score(np.asarray(y_true, dtype=np.int64), np.asarray(y_pred, dtype=np.int64))
    )


def _select_enet_hyperparams(
    x_train: np.ndarray,
    y_train: np.ndarray,
    *,
    task: Literal["age", "tissue", "sex"],
    study_ids: np.ndarray | None,
    alpha_grid: tuple[float, ...] = NESTED_ENET_ALPHA_GRID,
    l1_grid: tuple[float, ...] = NESTED_ENET_L1_GRID,
    seed: int = 42,
) -> tuple[float, float, dict[str, Any]]:
    """Pick α / l1_ratio on one study-grouped inner validation split."""
    n = int(x_train.shape[0])
    if study_ids is None:
        study_ids = np.asarray(["NA"] * n, dtype=object)
    splits = _study_inner_folds(study_ids, n_inner_folds=3, seed=seed)
    # Use the first non-degenerate split (plan: one inner validation split).
    inner_tr, inner_va = splits[0]
    if inner_tr.size < 2 or inner_va.size < 1:
        return float(alpha_grid[0]), float(l1_grid[0]), {"fallback": "tiny_inner_split"}

    best_score = -float("inf")
    best_alpha = float(alpha_grid[0])
    best_l1 = float(l1_grid[0])
    trials: list[dict[str, Any]] = []
    x_inner_tr, x_inner_va, _ = _scale_train_test(x_train[inner_tr], x_train[inner_va])
    y_tr = np.asarray(y_train)[inner_tr]
    y_va = np.asarray(y_train)[inner_va]
    sk_task: Literal["regression", "multiclass"] = "regression" if task == "age" else "multiclass"
    for alpha in alpha_grid:
        for l1_ratio in l1_grid:
            if sk_task == "multiclass" and len(np.unique(y_tr)) < 2:
                continue
            model = _fit_nested_enet_model(
                x_inner_tr,
                y_tr,
                task=sk_task,
                alpha=float(alpha),
                l1_ratio=float(l1_ratio),
            )
            if task == "age":
                pred = np.asarray(model.predict(x_inner_va), dtype=np.float64)
                score = _inner_val_score(task=task, y_true=y_va, y_pred=pred)
                proba = None
            else:
                pred = np.asarray(model.predict(x_inner_va), dtype=np.int64)
                proba = None
                if hasattr(model, "predict_proba"):
                    try:
                        proba = np.asarray(model.predict_proba(x_inner_va), dtype=np.float64)
                    except Exception:
                        proba = None
                score = _inner_val_score(task=task, y_true=y_va, y_pred=pred, y_proba=proba)
            trials.append(
                {"alpha": float(alpha), "l1_ratio": float(l1_ratio), "score": float(score)}
            )
            if score > best_score:
                best_score = score
                best_alpha = float(alpha)
                best_l1 = float(l1_ratio)
    meta = {
        "n_inner_train": int(inner_tr.size),
        "n_inner_val": int(inner_va.size),
        "best_score": float(best_score) if best_score > -float("inf") else None,
        "n_trials": len(trials),
    }
    return best_alpha, best_l1, meta


def run_nested_elasticnet_multitask(
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
    study_ids_train: np.ndarray | None = None,
    study_ids_test: np.ndarray | None = None,
    platforms_test: np.ndarray | None = None,
    tissue_class_names: list[str] | None = None,
    alpha_grid: tuple[float, ...] = NESTED_ENET_ALPHA_GRID,
    l1_grid: tuple[float, ...] = NESTED_ENET_L1_GRID,
    seed: int = 42,
) -> dict[str, Any]:
    """Nested elastic-net: train-fold StandardScaler + inner-val α/l1 selection.

    Unlike ``run_elasticnet_multitask`` (fixed α=0.1 / l1=0.5, no scaling), this
    selects hyperparameters on a study-grouped inner split of the outer train fold,
    then refits on the full outer train. Scaler is never fit on test rows.
    """
    x_tr_raw = np.asarray(x_train, dtype=np.float64)
    x_te_raw = np.asarray(x_test, dtype=np.float64)
    alpha_grid, l1_grid = _nested_grids_for_width(
        int(x_tr_raw.shape[1]), alpha_grid=alpha_grid, l1_grid=l1_grid
    )
    x_tr, x_te, scaler = _scale_train_test(x_tr_raw, x_te_raw)
    # Record train-only scaler stats for leak tests / diagnostics.
    scaler_mean = np.asarray(scaler.mean_, dtype=np.float64)
    models: dict[str, Any] = {}
    selected: dict[str, Any] = {}

    if (
        age_train is not None
        and age_mask_train is not None
        and bool(np.asarray(age_mask_train).any())
    ):
        m = np.asarray(age_mask_train, dtype=bool)
        studies = None if study_ids_train is None else np.asarray(study_ids_train)[m]
        alpha, l1_ratio, meta = _select_enet_hyperparams(
            x_tr_raw[m],
            np.asarray(age_train)[m],
            task="age",
            study_ids=studies,
            alpha_grid=alpha_grid,
            l1_grid=l1_grid,
            seed=seed,
        )
        models["age"] = _fit_nested_enet_model(
            x_tr[m],
            np.asarray(age_train)[m],
            task="regression",
            alpha=alpha,
            l1_ratio=l1_ratio,
        )
        selected["age"] = {"alpha": alpha, "l1_ratio": l1_ratio, **meta}

    if (
        tissue_train is not None
        and tissue_mask_train is not None
        and bool(np.asarray(tissue_mask_train).any())
    ):
        m = np.asarray(tissue_mask_train, dtype=bool)
        studies = None if study_ids_train is None else np.asarray(study_ids_train)[m]
        alpha, l1_ratio, meta = _select_enet_hyperparams(
            x_tr_raw[m],
            np.asarray(tissue_train)[m],
            task="tissue",
            study_ids=studies,
            alpha_grid=alpha_grid,
            l1_grid=l1_grid,
            seed=seed + 1,
        )
        models["tissue"] = _fit_nested_enet_model(
            x_tr[m],
            np.asarray(tissue_train)[m],
            task="multiclass",
            alpha=alpha,
            l1_ratio=l1_ratio,
        )
        selected["tissue"] = {"alpha": alpha, "l1_ratio": l1_ratio, **meta}

    if (
        sex_train is not None
        and sex_mask_train is not None
        and bool(np.asarray(sex_mask_train).any())
    ):
        m = np.asarray(sex_mask_train, dtype=bool)
        studies = None if study_ids_train is None else np.asarray(study_ids_train)[m]
        alpha, l1_ratio, meta = _select_enet_hyperparams(
            x_tr_raw[m],
            np.asarray(sex_train)[m],
            task="sex",
            study_ids=studies,
            alpha_grid=alpha_grid,
            l1_grid=l1_grid,
            seed=seed + 2,
        )
        models["sex"] = _fit_nested_enet_model(
            x_tr[m],
            np.asarray(sex_train)[m],
            task="multiclass",
            alpha=alpha,
            l1_ratio=l1_ratio,
        )
        selected["sex"] = {"alpha": alpha, "l1_ratio": l1_ratio, **meta}

    preds = predict_linear_multitask(models, x_te)
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
        "kind": "elasticnet_nested",
        "metrics": metrics,
        "n_features": int(x_tr_raw.shape[1]),
        "selected": selected,
        "scaler": {
            "with_mean": True,
            "fit_on": "outer_train_only",
            "n_train": int(x_tr_raw.shape[0]),
            "mean_l2": float(np.linalg.norm(scaler_mean)),
            "train_feature_mean": scaler_mean.tolist(),
        },
        "alpha_grid": list(alpha_grid),
        "l1_ratio_grid": list(l1_grid),
    }
