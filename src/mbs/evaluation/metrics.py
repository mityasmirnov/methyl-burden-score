"""Evaluation metrics for Stage 0 phenotype tasks (Milestone 5b)."""

from __future__ import annotations

from typing import Any

import numpy as np


def regression_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    """Age-style regression: MAE and RMSE (plus Pearson when defined)."""
    yt = np.asarray(y_true, dtype=np.float64).reshape(-1)
    yp = np.asarray(y_pred, dtype=np.float64).reshape(-1)
    if yt.shape != yp.shape:
        raise ValueError("y_true and y_pred shape mismatch")
    if yt.size == 0:
        raise ValueError("empty arrays")
    err = yp - yt
    mae = float(np.mean(np.abs(err)))
    rmse = float(np.sqrt(np.mean(err**2)))
    out: dict[str, float] = {"mae": mae, "rmse": rmse}
    if yt.size >= 2 and float(np.std(yt)) > 0:
        ss_tot = float(np.sum((yt - yt.mean()) ** 2))
        ss_res = float(np.sum(err**2))
        out["r2"] = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
        if float(np.std(yp)) > 0:
            out["pearson_r"] = float(np.corrcoef(yt, yp)[0, 1])
            ranks_t = np.argsort(np.argsort(yt))
            ranks_p = np.argsort(np.argsort(yp))
            out["spearman_r"] = float(np.corrcoef(ranks_t, ranks_p)[0, 1])
    return out


def _binary_ranks(y_true: np.ndarray, y_score: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    yt = np.asarray(y_true).reshape(-1)
    ys = np.asarray(y_score, dtype=np.float64).reshape(-1)
    if yt.shape != ys.shape:
        raise ValueError("y_true and y_score shape mismatch")
    classes = np.unique(yt)
    if classes.size != 2:
        raise ValueError(f"binary metrics require exactly 2 classes, got {classes.tolist()}")
    if set(classes.tolist()) != {0, 1}:
        pos = classes.max()
        yt = (yt == pos).astype(np.int64)
    else:
        yt = yt.astype(np.int64)
    return yt, ys


def binary_auroc_auprc(y_true: np.ndarray, y_score: np.ndarray) -> dict[str, float]:
    """AUROC and AUPRC (sklearn when available; numpy fallback otherwise)."""
    yt, ys = _binary_ranks(y_true, y_score)
    try:
        from sklearn.metrics import average_precision_score, roc_auc_score  # noqa: PLC0415

        return {
            "auroc": float(roc_auc_score(yt, ys)),
            "auprc": float(average_precision_score(yt, ys)),
        }
    except ImportError:
        # Mann-Whitney AUROC
        pos = ys[yt == 1]
        neg = ys[yt == 0]
        if pos.size == 0 or neg.size == 0:
            raise ValueError("need both positive and negative labels") from None
        # ponytail: O(n^2) pairwise AUROC; fine for Stage 0 unit sizes
        correct = 0.0
        for p in pos:
            correct += float(np.sum(p > neg)) + 0.5 * float(np.sum(p == neg))
        auroc = correct / (pos.size * neg.size)
        order = np.argsort(-ys)
        yt_sorted = yt[order]
        tp = 0
        fp = 0
        n_pos = int(yt.sum())
        precisions: list[float] = []
        recalls: list[float] = []
        for label in yt_sorted:
            if label == 1:
                tp += 1
            else:
                fp += 1
            precisions.append(tp / (tp + fp))
            recalls.append(tp / n_pos if n_pos else 0.0)
        if len(recalls) > 1:
            auprc = float(np.trapezoid(precisions, recalls))
        else:
            auprc = float(precisions[0])
        return {"auroc": float(auroc), "auprc": auprc}


def expected_calibration_error(
    y_true: np.ndarray,
    y_prob: np.ndarray,
    *,
    n_bins: int = 10,
) -> dict[str, float]:
    """Binary ECE over equal-width probability bins."""
    yt = np.asarray(y_true).reshape(-1).astype(np.float64)
    yp = np.asarray(y_prob).reshape(-1).astype(np.float64)
    if yt.shape != yp.shape or yt.size == 0:
        raise ValueError("y_true and y_prob length mismatch")
    bins = np.linspace(0.0, 1.0, n_bins + 1)
    ece = 0.0
    for i in range(n_bins):
        lo, hi = bins[i], bins[i + 1]
        mask = (yp >= lo) & (yp < hi) if i < n_bins - 1 else (yp >= lo) & (yp <= hi)
        if not mask.any():
            continue
        acc = float(yt[mask].mean())
        conf = float(yp[mask].mean())
        ece += (float(mask.mean())) * abs(acc - conf)
    return {"ece": float(ece)}


def multiclass_metrics(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    *,
    n_classes: int | None = None,
) -> dict[str, Any]:
    """Macro-F1, balanced accuracy, and confusion matrix (list-of-lists)."""
    yt = np.asarray(y_true).reshape(-1).astype(np.int64)
    yp = np.asarray(y_pred).reshape(-1).astype(np.int64)
    if yt.shape != yp.shape:
        raise ValueError("y_true and y_pred shape mismatch")
    if yt.size == 0:
        raise ValueError("empty arrays")
    k = int(n_classes) if n_classes is not None else int(max(yt.max(), yp.max()) + 1)
    confusion = np.zeros((k, k), dtype=np.int64)
    for t, p in zip(yt, yp, strict=True):
        if 0 <= t < k and 0 <= p < k:
            confusion[t, p] += 1
    recalls = []
    precisions = []
    f1s = []
    for c in range(k):
        tp = float(confusion[c, c])
        fn = float(confusion[c, :].sum() - tp)
        fp = float(confusion[:, c].sum() - tp)
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
        recalls.append(recall)
        precisions.append(precision)
        f1s.append(f1)
    return {
        "macro_f1": float(np.mean(f1s)),
        "balanced_accuracy": float(np.mean(recalls)),
        "confusion_matrix": confusion.tolist(),
        "per_class_precision": precisions,
        "per_class_recall": recalls,
        "per_class_f1": f1s,
    }


def metrics_by_group(
    y_true: np.ndarray,
    y_pred_or_score: np.ndarray,
    groups: np.ndarray,
    *,
    task: str,
    n_classes: int | None = None,
) -> dict[str, dict[str, Any]]:
    """Compute metrics stratified by study/platform group label."""
    groups = np.asarray(groups).reshape(-1)
    yt = np.asarray(y_true).reshape(-1)
    yp = np.asarray(y_pred_or_score).reshape(-1)
    if groups.shape[0] != yt.shape[0]:
        raise ValueError("groups length mismatch")
    out: dict[str, dict[str, Any]] = {}
    for g in sorted(set(groups.tolist()), key=str):
        mask = groups == g
        if task == "regression":
            out[str(g)] = regression_metrics(yt[mask], yp[mask])
        elif task == "binary":
            out[str(g)] = binary_auroc_auprc(yt[mask], yp[mask])
        elif task == "multiclass":
            out[str(g)] = multiclass_metrics(yt[mask], yp[mask], n_classes=n_classes)
        else:
            raise ValueError(f"unknown task: {task}")
    return out
