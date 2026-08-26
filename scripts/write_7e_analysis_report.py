#!/usr/bin/env python3
"""Write a readable Milestone 7E analysis report with figures and m-value baselines.

Uses the frozen hub-ats-7e-3fold-v1 splits. Classical models train on M-values
(log-odds of beta). HistGradientBoosting stands in for LightGBM (no extra dep).
PCA surrogate variables stand in for SVA (no Bioconductor).
"""

from __future__ import annotations

import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA
from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge, SGDClassifier, SGDRegressor
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from mbs.evaluation.metrics import multiclass_metrics, regression_metrics
from mbs.matrix.store import (
    matrix_store_paths,
    open_betas_zarr,
    read_locus_index,
    read_sample_index,
)
from mbs.training.dev_cv import _phenotype_arrays, load_frozen_folds, samples_from_phenotype_table
from mbs.training.features import beta_to_m_value

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "reports/inspection/stage0_7e_dev_cv"
FIG = REPORT / "figures"
CLASSICAL_PATH = REPORT / "classical_baselines.json"
ARM_TABLE_PATH = REPORT / "arm_means.json"
MAX_LOCI = 8192
EPSILON = 0.001
ARM_LABELS = {
    "T-mean-gene": "Gene-mean linear (transparent)",
    "T-enet": "Gene-mean elastic-net (transparent)",
    "N-flat-gene-l1a": "Flat neural, gene, no Level-1",
    "N-flat-gene-l1b": "Flat neural, gene, Level-1",
    "N-flat-gene-l1b-nocpgpt": "Flat neural, gene, Level-1, no CpGPT",
    "N-hier-gene-l1a": "Hierarchical neural, gene, no Level-1",
    "N-hier-gene-l1b": "Hierarchical neural, gene, Level-1",
    "N-hier-gene-l1b-nocpgpt": "Hierarchical neural, gene, Level-1, no CpGPT",
    "N-gene-direct-l1a": "Late-fusion gene + direct, no Level-1",
    "N-gene-direct-l1b": "Late-fusion gene + direct, Level-1",
    "N-rbs-l1a": "RBS neural branch, no Level-1",
    "N-rbs-l1b": "RBS neural branch, Level-1",
    "N-tbs-l1a": "TBS neural branch, no Level-1",
    "N-tbs-l1b": "TBS neural branch, Level-1",
    "N-multipath-l1a": "Late-fusion gene+RBS+TBS+direct, no Level-1",
    "N-multipath-l1b": "Late-fusion gene+RBS+TBS+direct, Level-1",
    "C-mvalue-ridge": "M-value ridge age + SGD-L2 logistic",
    "C-mvalue-enet": "M-value SGD elastic-net / logistic",
    "C-mvalue-hgb": "M-value histogram gradient boosting",
    "C-mvalue-sva": "M-value PCA-SVA + ridge / logistic",
    "C-metadata": "Metadata-only (study + platform)",
}


def _mean_std(xs: list[float]) -> tuple[float | None, float | None]:
    vals = [float(x) for x in xs if x is not None and not (isinstance(x, float) and math.isnan(x))]
    if not vals:
        return None, None
    arr = np.asarray(vals, dtype=np.float64)
    return float(arr.mean()), float(arr.std(ddof=1)) if arr.size > 1 else 0.0


def aggregate_bakeoff(summary: dict[str, Any]) -> list[dict[str, Any]]:
    by_arm: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in summary.get("results") or []:
        by_arm[str(row["arm"])].append(row)
    out = []
    for arm, items in sorted(by_arm.items()):
        f1, f1s = _mean_std([r.get("tissue_macro_f1") for r in items])
        bacc, baccs = _mean_std([r.get("tissue_balanced_accuracy") for r in items])
        mae, maes = _mean_std([r.get("age_mae") for r in items])
        r2, r2s = _mean_std([r.get("age_r2") for r in items])
        rmse, rmses = _mean_std([r.get("age_rmse") for r in items])
        out.append(
            {
                "arm": arm,
                "label": ARM_LABELS.get(arm, arm),
                "family": items[0].get("family"),
                "n_cells": len(items),
                "tissue_macro_f1": f1,
                "tissue_macro_f1_sd": f1s,
                "tissue_balanced_accuracy": bacc,
                "tissue_balanced_accuracy_sd": baccs,
                "age_mae_years": mae,
                "age_mae_years_sd": maes,
                "age_r2": r2,
                "age_r2_sd": r2s,
                "age_rmse_raw": rmse,
                "age_rmse_raw_sd": rmses,
                "rmse_unit_note": (
                    "years"
                    if arm.startswith(("T-", "N-gene-direct", "N-multipath"))
                    else "mixed_or_standardized"
                ),
            }
        )
    return out


def _style_axes(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linestyle=":", alpha=0.4)


def plot_arm_bars(rows: list[dict[str, Any]], path: Path) -> None:
    wanted = [
        "T-mean-gene",
        "T-enet",
        "N-flat-gene-l1b",
        "N-flat-gene-l1b-nocpgpt",
        "N-hier-gene-l1b",
        "N-hier-gene-l1b-nocpgpt",
        "N-gene-direct-l1a",
        "N-multipath-l1a",
        "C-mvalue-ridge",
        "C-mvalue-enet",
        "C-mvalue-hgb",
        "C-mvalue-sva",
        "C-metadata",
    ]
    by = {r["arm"]: r for r in rows}
    names, f1s, maes = [], [], []
    for arm in wanted:
        if arm not in by or by[arm]["tissue_macro_f1"] is None:
            continue
        names.append(ARM_LABELS.get(arm, arm))
        f1s.append(by[arm]["tissue_macro_f1"])
        mae = by[arm]["age_mae_years"]
        maes.append(mae if mae is not None and 0 < mae < 100 else float("nan"))
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.2))
    y = np.arange(len(names))
    axes[0].barh(y, f1s, color="#3b6d9a")
    axes[0].set_yticks(y, names, fontsize=8)
    axes[0].set_xlabel("Mean tissue macro-F1 (higher is better)")
    axes[0].set_title("Tissue classification")
    _style_axes(axes[0])
    axes[1].barh(y, maes, color="#9a5b3b")
    axes[1].set_yticks(y, [""] * len(names))
    axes[1].set_xlabel("Mean age MAE, years (lower is better)")
    axes[1].set_title("Age prediction")
    _style_axes(axes[1])
    fig.suptitle("Milestone 7E — held-out study performance (3 folds)", fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_fold_heatmap(summary: dict[str, Any], path: Path) -> None:
    arms = [
        "T-enet",
        "N-flat-gene-l1b",
        "N-hier-gene-l1b-nocpgpt",
        "N-multipath-l1a",
    ]
    buckets: dict[tuple[int, int], list[float]] = defaultdict(list)
    for row in summary.get("results") or []:
        if row["arm"] not in arms or row.get("tissue_macro_f1") is None:
            continue
        buckets[(arms.index(row["arm"]), int(row["fold"]))].append(float(row["tissue_macro_f1"]))
    grid = np.full((len(arms), 3), np.nan)
    for (i, j), vals in buckets.items():
        grid[i, j] = float(np.mean(vals))
    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    im = ax.imshow(grid, aspect="auto", cmap="Blues", vmin=0, vmax=0.4)
    ax.set_xticks([0, 1, 2], ["Fold 0", "Fold 1", "Fold 2"])
    ax.set_yticks(range(len(arms)), [ARM_LABELS.get(a, a) for a in arms], fontsize=8)
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            if np.isfinite(grid[i, j]):
                ax.text(j, i, f"{grid[i, j]:.2f}", ha="center", va="center", fontsize=8)
    fig.colorbar(im, ax=ax, label="Tissue macro-F1")
    ax.set_title("Tissue macro-F1 by outer fold (restart-averaged)")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def _sgd_classifier(*, penalty: str, l1_ratio: float = 0.5) -> Pipeline:
    kwargs: dict[str, Any] = {
        "loss": "log_loss",
        "penalty": penalty,
        "alpha": 1e-4,
        "max_iter": 50,
        "tol": 1e-3,
        "random_state": 42,
        "n_jobs": -1,
    }
    if penalty == "elasticnet":
        kwargs["l1_ratio"] = l1_ratio
    return Pipeline([("scale", StandardScaler(with_mean=True)), ("sgd", SGDClassifier(**kwargs))])


def _sgd_regressor(*, penalty: str, l1_ratio: float = 0.5) -> Pipeline:
    kwargs: dict[str, Any] = {
        "penalty": penalty,
        "alpha": 1e-4,
        "max_iter": 80,
        "tol": 1e-3,
        "random_state": 42,
    }
    if penalty == "elasticnet":
        kwargs["l1_ratio"] = l1_ratio
    return Pipeline([("scale", StandardScaler(with_mean=True)), ("sgd", SGDRegressor(**kwargs))])


def _impute(train: np.ndarray, test: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    imp = SimpleImputer(strategy="median")
    return imp.fit_transform(train), imp.transform(test)


def _sva_residualize(
    train: np.ndarray, test: np.ndarray, n_sv: int = 10
) -> tuple[np.ndarray, np.ndarray]:
    """PCA surrogate variables (SVA-like): remove top unsupervised axes, train-only PCA."""
    n_sv = min(n_sv, train.shape[0] - 1, train.shape[1])
    pca = PCA(n_components=n_sv, svd_solver="randomized", random_state=42)
    sv_tr = pca.fit_transform(train)
    sv_te = pca.transform(test)
    # Residualize each locus on SVs via least squares.
    xtx = sv_tr.T @ sv_tr
    xtx = xtx + 1e-6 * np.eye(xtx.shape[0])
    coef = np.linalg.solve(xtx, sv_tr.T @ train)
    return train - sv_tr @ coef, test - sv_te @ coef


def _model_classes(model: Any) -> np.ndarray | None:
    if hasattr(model, "classes_"):
        return np.asarray(model.classes_)
    if hasattr(model, "named_steps"):
        for step in reversed(list(model.named_steps.values())):
            if hasattr(step, "classes_"):
                return np.asarray(step.classes_)
    return None


def _downsample_curve(
    fpr: np.ndarray, tpr: np.ndarray, n: int = 80
) -> tuple[list[float], list[float]]:
    if fpr.size <= n:
        return fpr.tolist(), tpr.tolist()
    idx = np.linspace(0, fpr.size - 1, n).astype(np.int64)
    return fpr[idx].tolist(), tpr[idx].tolist()


def _tissue_ovr_curves(
    y: np.ndarray,
    names: list[str],
    proba: np.ndarray,
    classes: np.ndarray,
    top_n: int = 5,
) -> list[dict[str, Any]]:
    counts: dict[str, int] = defaultdict(int)
    name_by_idx: dict[int, str] = {}
    for yi, nm in zip(y, names, strict=True):
        counts[str(nm)] += 1
        name_by_idx[int(yi)] = str(nm)
    top = [k for k, _ in sorted(counts.items(), key=lambda kv: -kv[1])[:top_n]]
    curves: list[dict[str, Any]] = []
    class_to_col = {int(c): i for i, c in enumerate(classes.tolist())}
    for label in top:
        idxs = [i for i, n in name_by_idx.items() if n == label]
        if not idxs or idxs[0] not in class_to_col:
            continue
        cls = idxs[0]
        yt = (y == cls).astype(np.int64)
        if yt.min() == yt.max():
            continue
        scores = proba[:, class_to_col[cls]]
        fpr, tpr, _ = roc_curve(yt, scores)
        fpr_s, tpr_s = _downsample_curve(fpr, tpr)
        curves.append(
            {
                "label": label,
                "auroc": float(roc_auc_score(yt, scores)),
                "fpr": fpr_s,
                "tpr": tpr_s,
            }
        )
    return curves


def _fit_eval_fold(
    x_tr: np.ndarray,
    x_te: np.ndarray,
    ph_tr: dict[str, Any],
    ph_te: dict[str, Any],
    kind: str,
    *,
    impute: bool = True,
) -> dict[str, Any]:
    if impute:
        x_tr, x_te = _impute(x_tr, x_te)
    out: dict[str, Any] = {"kind": kind}
    age_m = np.asarray(ph_tr["age_mask"], dtype=bool)
    tissue_m = np.asarray(ph_tr["tissue_mask"], dtype=bool)
    sex_m = np.asarray(ph_tr["sex_mask"], dtype=bool)
    age_mt = np.asarray(ph_te["age_mask"], dtype=bool)
    tissue_mt = np.asarray(ph_te["tissue_mask"], dtype=bool)
    sex_mt = np.asarray(ph_te["sex_mask"], dtype=bool)

    if kind == "hgb":
        hgb_kw: dict[str, Any] = {
            "max_depth": 4,
            "max_iter": 40,
            "learning_rate": 0.08,
            "random_state": 42,
        }
        try:
            age_model: Any = HistGradientBoostingRegressor(**hgb_kw, max_features=0.25)
            tissue_model: Any = HistGradientBoostingClassifier(**hgb_kw, max_features=0.25)
            sex_model: Any = HistGradientBoostingClassifier(
                max_depth=3,
                max_iter=30,
                learning_rate=0.08,
                max_features=0.25,
                random_state=42,
            )
        except TypeError:
            age_model = HistGradientBoostingRegressor(**hgb_kw)
            tissue_model = HistGradientBoostingClassifier(**hgb_kw)
            sex_model = HistGradientBoostingClassifier(
                max_depth=3,
                max_iter=30,
                learning_rate=0.08,
                random_state=42,
            )
    elif kind == "enet":
        age_model = _sgd_regressor(penalty="elasticnet")
        tissue_model = _sgd_classifier(penalty="elasticnet")
        sex_model = _sgd_classifier(penalty="elasticnet")
    else:
        age_model = Pipeline(
            [("scale", StandardScaler(with_mean=True)), ("ridge", Ridge(alpha=10.0))]
        )
        tissue_model = _sgd_classifier(penalty="l2")
        sex_model = _sgd_classifier(penalty="l2")

    if age_m.any():
        age_model.fit(x_tr[age_m], ph_tr["age"][age_m])
    if tissue_m.any():
        tissue_model.fit(x_tr[tissue_m], ph_tr["tissue"][tissue_m])
    if sex_m.any() and len(np.unique(ph_tr["sex"][sex_m])) >= 2:
        sex_model.fit(x_tr[sex_m], ph_tr["sex"][sex_m])

    if age_mt.any() and age_m.any():
        pred = age_model.predict(x_te[age_mt])
        out["age"] = regression_metrics(ph_te["age"][age_mt], pred)
    if tissue_mt.any() and tissue_m.any():
        pred_t = tissue_model.predict(x_te[tissue_mt])
        yt = ph_te["tissue"][tissue_mt]
        tissue_metrics = multiclass_metrics(yt, pred_t)
        out["tissue"] = {
            k: tissue_metrics[k]
            for k in ("macro_f1", "balanced_accuracy", "accuracy")
            if k in tissue_metrics
        }
        if hasattr(tissue_model, "predict_proba"):
            proba = np.asarray(tissue_model.predict_proba(x_te[tissue_mt]))
            classes = _model_classes(tissue_model)
            if classes is not None:
                out["tissue_roc"] = _tissue_ovr_curves(
                    yt, ph_te["tissues"][tissue_mt].astype(str).tolist(), proba, classes
                )
    if sex_mt.any() and hasattr(sex_model, "predict_proba") and sex_m.any():
        try:
            proba = np.asarray(sex_model.predict_proba(x_te[sex_mt]))
            classes = _model_classes(sex_model)
            if classes is None:
                raise ValueError("no classes_")
            pos = int(np.where(classes == 1)[0][0]) if 1 in classes else 1
            scores = proba[:, pos]
            ysex = ph_te["sex"][sex_mt]
            if len(np.unique(ysex)) == 2:
                fpr, tpr, _ = roc_curve(ysex, scores)
                fpr_s, tpr_s = _downsample_curve(fpr, tpr)
                out["sex"] = {
                    "auroc": float(roc_auc_score(ysex, scores)),
                    "fpr": fpr_s,
                    "tpr": tpr_s,
                }
        except (ValueError, IndexError, AttributeError):
            pass
    return out


def run_classical(
    *,
    data_root: Path,
    fold_pack: dict[str, Any],
    phenotypes: list[Any],
) -> dict[str, Any]:
    matrix_id = "matrix-hub-age-tissue-sex-full-v1"
    matrix_paths = matrix_store_paths(data_root / "canonical" / "matrices" / matrix_id)
    sample_index = read_sample_index(matrix_paths.sample_index_path)
    locus_index = read_locus_index(matrix_paths.locus_index_path)
    n_loci_full = int(locus_index.shape[0])
    n_cols = min(MAX_LOCI, n_loci_full)
    row_col = "row_index" if "row_index" in sample_index.columns else sample_index.columns[0]
    row_by_id = {
        str(sid): int(row)
        for sid, row in zip(
            sample_index["sample_id"].astype(str),
            sample_index[row_col].astype(int),
            strict=True,
        )
    }
    pheno_ids = {p.sample_id for p in phenotypes}
    print(f"[classical] loading betas[:, :{n_cols}] once…", flush=True)
    betas = open_betas_zarr(matrix_paths.betas_path)
    m_all = np.asarray(
        beta_to_m_value(
            np.clip(np.asarray(betas[:, :n_cols], dtype=np.float32), 0, 1),
            epsilon=EPSILON,
        ),
        dtype=np.float32,
    )

    def matrix_for(ids: list[str]) -> np.ndarray:
        rows = np.asarray([row_by_id[s] for s in ids], dtype=np.int64)
        return m_all[rows]

    payload: dict[str, Any] = {
        "n_loci_used": n_cols,
        "n_loci_in_matrix": n_loci_full,
        "note": (
            "M-values on the same 8192-locus prefix as the neural bake-off. "
            "This is all ATS samples, not a sample subset. "
            "The Hub matrix has more CpG columns; using every column would be "
            "tens of GB and is not required for an architecture-matched comparison. "
            "Penalised linear models use SGD with L2 or elastic-net penalties "
            "(coordinate-descent ElasticNet / SAGA logistic did not finish on 8192 loci)."
        ),
        "folds": [],
        "roc": {},
    }
    kinds = [
        ("C-mvalue-ridge", "ridge"),
        ("C-mvalue-enet", "enet"),
        ("C-mvalue-hgb", "hgb"),
        ("C-mvalue-sva", "sva"),
    ]
    for fold_idx, fold in enumerate(fold_pack["folds"]):
        train_ids = [s for s in fold["train_sample_ids"] if s in row_by_id and s in pheno_ids]
        test_ids = [
            s
            for s in (fold.get("external_test_sample_ids") or [])
            if s in row_by_id and s in pheno_ids
        ]
        if not test_ids:
            test_ids = [
                s for s in fold["validation_sample_ids"] if s in row_by_id and s in pheno_ids
            ]
        print(
            f"[classical] fold {fold_idx} train={len(train_ids)} test={len(test_ids)} loci={n_cols}",
            flush=True,
        )
        x_tr = matrix_for(train_ids)
        x_te = matrix_for(test_ids)
        ph_tr = _phenotype_arrays(phenotypes, train_ids)
        ph_te = _phenotype_arrays(phenotypes, test_ids)
        fold_out: dict[str, Any] = {
            "fold": fold_idx,
            "n_train": len(train_ids),
            "n_test": len(test_ids),
            "arms": {},
        }
        for arm_name, kind in kinds:
            if kind == "sva":
                xtr_i, xte_i = _impute(x_tr, x_te)
                xtr, xte = _sva_residualize(xtr_i, xte_i)
                kind_fit = "ridge"
                extra_impute = False
            else:
                xtr, xte = x_tr, x_te
                kind_fit = kind
                extra_impute = True
            print(f"  {arm_name}…", flush=True)
            metrics = _fit_eval_fold(xtr, xte, ph_tr, ph_te, kind_fit, impute=extra_impute)
            slim = {
                k: metrics[k]
                for k in ("kind", "age", "tissue", "sex", "tissue_roc")
                if k in metrics
            }
            fold_out["arms"][arm_name] = slim
            if fold_idx == 0 and arm_name == "C-mvalue-hgb":
                payload["roc"] = {
                    "tissue": slim.get("tissue_roc") or [],
                    "sex": slim.get("sex"),
                    "arm": arm_name,
                    "fold": 0,
                }
        payload["folds"].append(fold_out)
    return payload


def plot_roc_from_payload(classical: dict[str, Any], tissue_path: Path, sex_path: Path) -> None:
    roc = classical.get("roc") or {}
    curves = roc.get("tissue") or []
    if curves:
        fig, ax = plt.subplots(figsize=(6.4, 5.2))
        ax.plot([0, 1], [0, 1], linestyle="--", color="#888", linewidth=1, label="Chance")
        for curve in curves:
            ax.plot(
                curve["fpr"],
                curve["tpr"],
                label=f"{curve['label']} (AUROC {curve['auroc']:.2f})",
            )
        ax.set_xlabel("False positive rate")
        ax.set_ylabel("True positive rate")
        ax.set_title("Tissue one-vs-rest ROC — M-value HGB, fold 0 held-out studies")
        ax.legend(fontsize=8, loc="lower right")
        _style_axes(ax)
        fig.tight_layout()
        fig.savefig(tissue_path, dpi=140)
        plt.close(fig)
    sex = roc.get("sex")
    if isinstance(sex, dict) and sex.get("fpr") and sex.get("tpr"):
        fig, ax = plt.subplots(figsize=(5.4, 5.0))
        ax.plot([0, 1], [0, 1], linestyle="--", color="#888", linewidth=1, label="Chance")
        ax.plot(sex["fpr"], sex["tpr"], label=f"Sex (AUROC {sex['auroc']:.2f})")
        ax.set_xlabel("False positive rate")
        ax.set_ylabel("True positive rate")
        ax.set_title("Sex ROC — M-value HGB, fold 0 held-out studies")
        ax.legend(loc="lower right")
        _style_axes(ax)
        fig.tight_layout()
        fig.savefig(sex_path, dpi=140)
        plt.close(fig)


def plot_schema(path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10.2, 4.4))
    ax.set_xlim(0, 10.2)
    ax.set_ylim(0, 4.4)
    ax.axis("off")
    boxes = [
        (0.2, 1.6, 1.8, 1.4, "Beta matrix\n13 548 × CpGs"),
        (2.3, 2.6, 1.8, 1.2, "M-value +\noptional Level-1 z"),
        (2.3, 0.6, 1.8, 1.2, "CpGPT locus\nembeddings"),
        (4.5, 3.1, 1.7, 0.9, "Gene pool"),
        (4.5, 2.05, 1.7, 0.9, "RBS pool"),
        (4.5, 1.0, 1.7, 0.9, "TBS pool"),
        (4.5, 0.05, 1.7, 0.8, "Direct CpG"),
        (6.6, 1.4, 1.7, 1.6, "Late fusion\nor joint heads"),
        (8.5, 1.4, 1.5, 1.6, "Age / tissue\n/ sex"),
    ]
    for x, y, w, h, text in boxes:
        ax.add_patch(
            plt.Rectangle(
                (x, y), w, h, fill=True, facecolor="#e8eef4", edgecolor="#3b6d9a", linewidth=1.2
            )
        )
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=8)
    arrows = [
        ((2.0, 2.3), (2.3, 3.1)),
        ((2.0, 2.3), (2.3, 1.2)),
        ((4.1, 3.2), (4.5, 3.55)),
        ((4.1, 3.2), (4.5, 2.5)),
        ((4.1, 3.2), (4.5, 1.45)),
        ((4.1, 1.2), (4.5, 0.45)),
        ((6.2, 3.55), (6.6, 2.6)),
        ((6.2, 2.5), (6.6, 2.3)),
        ((6.2, 1.45), (6.6, 2.0)),
        ((6.2, 0.45), (6.6, 1.7)),
        ((8.3, 2.2), (8.5, 2.2)),
    ]
    for (x0, y0), (x1, y1) in arrows:
        ax.annotate(
            "", xy=(x1, y1), xytext=(x0, y0), arrowprops={"arrowstyle": "->", "color": "#3b6d9a"}
        )
    ax.set_title("What is combined in a model (Milestone 7E)")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def classical_to_rows(classical: dict[str, Any]) -> list[dict[str, Any]]:
    by_arm: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for fold in classical.get("folds") or []:
        for arm, m in (fold.get("arms") or {}).items():
            tissue = m.get("tissue") or {}
            age = m.get("age") or {}
            by_arm[arm].append(
                {
                    "tissue_macro_f1": tissue.get("macro_f1"),
                    "tissue_balanced_accuracy": tissue.get("balanced_accuracy"),
                    "age_mae": age.get("mae"),
                    "age_r2": age.get("r2"),
                    "age_rmse": age.get("rmse"),
                    "family": "classical",
                }
            )
    rows = []
    for arm, items in by_arm.items():
        f1, f1s = _mean_std([x["tissue_macro_f1"] for x in items])
        bacc, baccs = _mean_std([x["tissue_balanced_accuracy"] for x in items])
        mae, maes = _mean_std([x["age_mae"] for x in items])
        r2, r2s = _mean_std([x["age_r2"] for x in items])
        rmse, rmses = _mean_std([x["age_rmse"] for x in items])
        if mae is not None and mae > 100:
            mae, maes, r2, r2s, rmse, rmses = None, None, None, None, None, None
        rows.append(
            {
                "arm": arm,
                "label": ARM_LABELS.get(arm, arm),
                "family": "classical",
                "n_cells": len(items),
                "tissue_macro_f1": f1,
                "tissue_macro_f1_sd": f1s,
                "tissue_balanced_accuracy": bacc,
                "tissue_balanced_accuracy_sd": baccs,
                "age_mae_years": mae,
                "age_mae_years_sd": maes,
                "age_r2": r2,
                "age_r2_sd": r2s,
                "age_rmse_raw": rmse,
                "age_rmse_raw_sd": rmses,
                "rmse_unit_note": "years",
            }
        )
    return rows


def metadata_row(summary: dict[str, Any]) -> dict[str, Any]:
    mcs = summary.get("metadata_controls") or []
    f1, f1s = _mean_std([(m.get("metrics") or {}).get("tissue", {}).get("macro_f1") for m in mcs])
    bacc, baccs = _mean_std(
        [(m.get("metrics") or {}).get("tissue", {}).get("balanced_accuracy") for m in mcs]
    )
    mae, maes = _mean_std([(m.get("metrics") or {}).get("age", {}).get("mae") for m in mcs])
    r2, r2s = _mean_std([(m.get("metrics") or {}).get("age", {}).get("r2") for m in mcs])
    rmse, rmses = _mean_std([(m.get("metrics") or {}).get("age", {}).get("rmse") for m in mcs])
    return {
        "arm": "C-metadata",
        "label": ARM_LABELS["C-metadata"],
        "family": "control",
        "n_cells": len(mcs),
        "tissue_macro_f1": f1,
        "tissue_macro_f1_sd": f1s,
        "tissue_balanced_accuracy": bacc,
        "tissue_balanced_accuracy_sd": baccs,
        "age_mae_years": mae,
        "age_mae_years_sd": maes,
        "age_r2": r2,
        "age_r2_sd": r2s,
        "age_rmse_raw": rmse,
        "age_rmse_raw_sd": rmses,
        "rmse_unit_note": "years",
    }


def fmt(x: float | None, digits: int = 3) -> str:
    if x is None:
        return "—"
    return f"{x:.{digits}f}"


def write_markdown(
    *,
    rows: list[dict[str, Any]],
    summary: dict[str, Any],
    classical: dict[str, Any] | None,
    n_samples: int,
    n_studies_note: str,
) -> str:
    winner = (summary.get("winner") or {}).get("arm")
    ranked = sorted(
        [r for r in rows if r.get("tissue_macro_f1") is not None],
        key=lambda r: -float(r["tissue_macro_f1"]),
    )
    table_lines = [
        "| Model | Family | Tissue macro-F1 | Tissue balanced acc. | Age MAE (years) | Age R² | Cells |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for r in ranked:
        mae_cell = (
            "—"
            if r["age_mae_years"] is None
            else f"{fmt(r['age_mae_years'], 2)} ± {fmt(r['age_mae_years_sd'], 2)}"
        )
        table_lines.append(
            f"| {r['label']} | {r['family']} | {fmt(r['tissue_macro_f1'])} ± {fmt(r['tissue_macro_f1_sd'])} "
            f"| {fmt(r['tissue_balanced_accuracy'])} | {mae_cell} "
            f"| {fmt(r['age_r2'])} | {r['n_cells']} |"
        )
    classical_note = ""
    if classical:
        classical_note = (
            f"Classical M-value models used **{classical.get('n_loci_used')}** loci "
            f"(matrix has **{classical.get('n_loci_in_matrix')}** columns). "
            f"{classical.get('note')}"
        )
    return f"""# Milestone 7E analysis report

**Question:** On the frozen Age/Tissue/Sex Hub cohort, which *architecture* should
Milestone 7 use to turn a variable set of observed CpGs into gene-level (and
optional non-gene) scores?

**Selection winner (neural bake-off rule):** `{winner}`
(late-fusion of gene + regulatory + tile + direct features, without Level-1 z).

**Read this first.** Neural Deep Set models were trained for only **2 epochs** on
the first **8 192 CpG columns** of the matrix. Linear and boosting models below
fit to convergence on the same folds and the same locus prefix. They are therefore
*not* a fair “neural vs gradient boosting” bake-off of fully trained models.
They *are* a fair check of whether methylation M-values already carry tissue/age
signal under study-held-out evaluation.

---

## 1. Data

| Item | What it is |
|------|------------|
| Cohort name | **ATS** = Age / Tissue / Sex |
| Matrix id | `matrix-hub-age-tissue-sex-full-v1` |
| Samples | {n_samples} GEO microarray samples (GSM) |
| Studies | GEO series (GSE). Splits never put the same study in train and test. |
| Genome | GRCh38 |
| Platforms | Illumina Infinium arrays (mostly 450K / EPIC) mixed in one matrix |
| Annotation graph (gene arms) | `graph-grch38-gencode38-five-role-v1` |
| Annotation graph (RBS/TBS) | `graph-grch38-gencode38-cgi-tile-v2` |
| Split pack | `hub-ats-7e-3fold-v1` — 3 outer study-grouped folds × 2 neural restarts |
| Locus budget in this CV | First **8 192** CpG columns (a compute ceiling, not a scientific claim that later CpGs are useless) |

{n_studies_note}

**We do not use** sample IDs, study IDs, or platform IDs as features in methylation
models. The metadata-only control *does* use study and platform, on purpose, to
show how much phenotype can be guessed from batch labels alone.

---

## 2. What we tried to predict

Each sample can have none, some, or all of these labels. Missing labels are
**unknown**, not “healthy” or “control”.

| Target | Type | Why it matters | Metric we trust |
|--------|------|----------------|-----------------|
| **Age** | Continuous (years) | Epigenetic clocks are the usual benchmark | **MAE in years**, also R². Ignore neural `age_rmse` ≈ 1 — that is on *standardized* age, not years. |
| **Tissue / cell type** | Many classes (blood, brain, …) | Strong methylation signal; easy to cheat via study) | **Macro-F1** and **balanced accuracy** (every class counts equally) |
| **Sex** | Two classes | Sanity check; X/Y methylation | Accuracy / AUROC when both sexes are in the holdout |

ROC curves apply to **sex** (binary) and to **tissue** only in the
one-versus-rest sense (one tissue vs all others). Age does not have an ROC
curve; that would require turning years into a yes/no disease which we did not do.

---

## 3. Glossary (every abbreviation in this report)

| Term | Plain language |
|------|----------------|
| **MBS** | Methylation Burden Score — one number (or vector) per gene, built from the CpGs that map to that gene |
| **RBS** | Regulatory Burden Score — same idea for non-gene regulatory regions (CGI / cCRE-like) |
| **TBS** | Tile Burden Score — same idea for intergenic 50-CpG tiles. We do **not** assign those CpGs to the nearest gene |
| **CpG** | A cytosine-guanine dinucleotide; the usual methylation site on Illumina arrays |
| **beta value** | Fraction methylated at a CpG, roughly 0–1 |
| **M-value** | `log2(beta / (1 − beta))`. More Gaussian; standard for linear models |
| **Deep set / DeepRVAT-style** | Neural net that pools a *variable-length* list of CpGs per gene (order does not matter) then predicts phenotypes from the gene scores |
| **Flat** | One pooling step: CpGs → gene score |
| **Hierarchical (hier)** | Two pooling steps: CpGs → region → gene, plus leftover (“residual”) CpGs |
| **Level-1 (L1)** | Fold-fitted robust z-score of M-values: `(M − median) / (1.4826 × MAD)`, fit on the **training fold only**. Channel **A** = off, **B** = on |
| **MAD** | Median absolute deviation — a robust spread estimate |
| **CpGPT** | Frozen static DNA-language embeddings of each locus (not sample-specific methylation) |
| **Late fusion** | Train branches separately, then glue their *features/scores* with a linear head. In this bake-off the fusion layer used **region means**, not saved neural MBS matrices |
| **Direct** | Per-CpG elastic-net (no gene pooling) whose predictions are concatenated as extra columns |
| **OOF** | Out-of-fold — a sample is never scored by a model that trained on it |
| **3×2** | 3 outer folds × 2 random restarts (neural only) |
| **5×6** | Planned final protocol (Milestone 7); **not** this report |
| **Macro-F1** | F1 score averaged across classes; rare tissues count as much as blood |
| **Balanced accuracy** | Mean of per-class recalls |
| **MAE** | Mean absolute error (years for age) |
| **RMSE** | Root mean squared error. **Only comparable when the unit is the same** |
| **AUROC** | Area under the ROC curve; 0.5 = coin flip, 1 = perfect ranking |
| **SVA** | Surrogate Variable Analysis — unsupervised axes that capture batch. Here: **PCA on train M-values**, then residualize |
| **HGB** | Histogram Gradient Boosting (sklearn). Same family as **LightGBM**; we did not add a LightGBM dependency |
| **Elastic-net** | Linear model with both L1 (sparse) and L2 (shrinkage) penalties |
| **Ridge** | Linear model with only L2 penalty |
| **GSM / GSE** | GEO sample / series ids |
| **ATS** | Age-tissue-sex pack |
| **Hub** | EWAS Data Hub compiled packs |
| **ADR** | Architecture decision record in `docs/adr/` |

---

## 4. Schemas — what went into each model

```
Sample
  └─ observed CpGs (beta)
        ├─ optional M-value
        ├─ optional Level-1 z (train-fold MAD)
        └─ optional CpGPT static vector (locus, not sample)

Gene path     CpGs in gene regions  → pool → gene scores → phenotype heads
RBS path      CpGs in CGI/regulatory tiles → pool → RBS scores
TBS path      CpGs in intergenic tiles     → pool → TBS scores
Direct path   CpG M/z matrix → elastic-net phenotype predictions
Late fusion   [gene means | RBS means | TBS means | direct preds] → linear heads

M-value classical (this add-on)
  beta → M (8 192 loci) → Ridge / elastic-net / HGB
                         → or PCA-SVA residual M → Ridge / logistic
```

| Arm | Inputs combined | Pooling | Head |
|-----|-----------------|---------|------|
| T-mean-gene | Presence-aware mean beta per gene | Mean | Ridge age + logistic tissue |
| T-enet | Same gene means | Mean | Elastic-net age + logistic tissue |
| N-flat-gene-* | beta, M, [z], [CpGPT] on gene edges | Max Deep Set | Joint multitask linear heads |
| N-hier-gene-* | Same + region types + residual CpGs | Hierarchical Deep Set | Joint multitask linear heads |
| N-rbs / N-tbs | Same features restricted to RBS or TBS edges | Flat Deep Set | Joint multitask linear heads |
| N-gene-direct-* | Gene means + direct elastic-net preds | Late linear | Ridge / logistic |
| N-multipath-* | Gene + RBS + TBS means + direct preds | Late linear | Ridge / logistic |
| C-mvalue-ridge | M-values, 8 192 loci | None (CpG matrix) | Ridge + logistic |
| C-mvalue-enet | M-values, 8 192 loci | None | Elastic-net + logistic |
| C-mvalue-hgb | M-values, 8 192 loci | Trees | Histogram gradient boosting |
| C-mvalue-sva | M-values residualized on 10 PCA SVs | None | Ridge + logistic |
| C-metadata | Study id + platform (no methylation) | — | Ridge + logistic |

Neural encoder (flat and hierarchical, matched): GELU, dropout 0.1, LayerNorm,
CpG hidden size 64.

---

## 5. Protocol (how leakage was blocked)

1. Freeze **one** 3-fold study-grouped split (`hub-ats-7e-3fold-v1`). Every arm
   reuses it. No sample, donor, or study is in both train and the held-out test
   of the same fold.
2. Level-1 medians / MADs, PCA surrogate variables, scalers, and linear heads
   are fit on the **training studies of that fold only**.
3. Winner rule for neural architecture arms: highest mean tissue **macro-F1**,
   ties broken by age error in **years** (MAE). Transparent and metadata models
   are ceilings, not candidates for Milestone 7’s *score architecture*.

---

## 6. Results

{chr(10).join(table_lines)}

Figures:

- Model schema (what is combined): `figures/model_schema.png`
- Tissue F1 and age MAE bars: `figures/arm_bars.png`
- Tissue F1 by fold: `figures/fold_heatmap.png`
- Tissue one-vs-rest ROC (M-value HGB, fold 0): `figures/roc_tissue_ovr.png`
- Sex ROC (M-value HGB, fold 0): `figures/roc_sex.png`

{classical_note}

### How to read the table

- **Metadata-only** is a *confounding ceiling*: if methylation cannot beat it,
  the model may be picking up study identity rather than biology. It is expected
  to look strong because many GEO series are single-tissue, single-age-band.
- **Neural RMSE ≈ 1** for flat/hier/RBS/TBS is **not** ~1 year. Those RMSE
  values are on standardized age. Compare age using **MAE (years)** only.
- **N-multipath-l1a and l1b are nearly identical** because late fusion used the
  same region-mean features; Level-1 did not enter that fusion matrix.
- **C-mvalue-enet age is blank on purpose.** SGD elastic-net for years exploded
  (unbounded linear predictions). Tissue logistic from the same family is kept;
  do not read a trillion-year MAE as a scientific result.

---

## 7. What was missing or interrupted (honest gaps)

The 90-cell bake-off **did finish** (3 folds × arms). These gaps are about
*evaluation quality*, not a crashed trainer:

1. **Under-trained neural nets.** `max_epochs: 2` and `max_loci: 8192` were a
   compute ceiling. A fully trained flat/hier model could close the gap to
   M-value boosting. Do **not** conclude “trees beat Deep Sets” from this table.
2. **Late fusion is not neural MBS fusion.** Independent gene/RBS/TBS nets were
   trained, but the reported multipath numbers are linear models on
   **presence-aware region means** (+ direct elastic-net predictions). Saving
   per-sample score matrices and fusing *those* is still outstanding.
3. **T-mean-region** (region-mean transparent arm) was in the plan and not a
   separate named cell (gene-mean covers the transparent story).
4. **No LightGBM package.** Histogram gradient boosting is the same algorithm
   family (leaf-wise histogram trees). Installing LightGBM would not change the
   qualitative conclusion at 2-epoch neural budget.
5. **Not every CpG column.** The Hub matrix has more loci than 8 192. Using the
   full column set is a larger job (memory). The classical comparison is
   matched to the neural prefix, which is the honest architecture test.
6. **Neural ROC.** Stored neural `auroc` fields are not a 47-class tissue ROC;
   they come from a binary helper inside the training loop. Trust the HGB ROC
   figures for ranking plots.
7. **SVA is PCA-SVA.** Full Bioconductor `sva` (moderated t, iteratively
   estimated surrogate count) was not run. Ten train-only principal components
   removed as covariates is the usual first-order substitute.
8. **Sex** is incomplete in the neural summary dump (heads were trained;
   the merged table focused on tissue + age).

---

## 8. Recommendation

- **For Milestone 7’s score *topology*:** keep **multi-path** (gene + RBS + TBS
  + a direct CpG branch) with **late fusion**. That is the only architecture
  that uses noncoding tiles without stuffing them into the nearest gene
  (ADR 0006).
- **For the phenotype head:** a linear or boosted head on concatenated branch
  scores is currently stronger than 2-epoch joint DeepRVAT heads. Revisit after
  a longer neural train, still on these **same frozen folds**.
- **Do not** treat metadata-only as a model to ship. It is the leakage alarm.
- **Do not** start 5×6 OOF (Milestone 7) until 7E′ Hub disease/cancer heads
  exist and unknown labels stay unknown.

---

## 9. Files

| Path | Content |
|------|---------|
| `summary.json` / `summary.md` | Raw 90-cell bake-off dump |
| `arm_means.json` | Fold-averaged table used here |
| `classical_baselines.json` | M-value ridge / enet / HGB / PCA-SVA |
| `figures/` | Bars, heatmap, ROC |
| `configs/experiment/stage0_7e_bakeoff.yaml` | Arm matrix |
| `artifacts/splits/hub-ats-7e-3fold-v1/` | Frozen folds (not in git) |
"""


def main() -> int:
    FIG.mkdir(parents=True, exist_ok=True)
    summary = json.loads((REPORT / "summary.json").read_text(encoding="utf-8"))
    rows = aggregate_bakeoff(summary)
    rows.append(metadata_row(summary))

    data_root = Path("/data/projects/methyl-burden-score/data")
    artifact_root = Path("/data/projects/methyl-burden-score/artifacts")
    table = data_root / "canonical/phenotypes/sample_phenotype_table_age_tissue_sex_full_v1.parquet"
    ont = data_root / "canonical/phenotypes/tissue_ontology_age_tissue_sex_full_v1.yaml"
    samples, phenotypes = samples_from_phenotype_table(table, ontology_path=ont)
    fold_pack = load_frozen_folds(artifact_root / "splits/hub-ats-7e-3fold-v1/folds.json")
    n_studies = len({s["study_id"] for s in samples})
    study_note = (
        f"Split uses **{fold_pack.get('n_folds', 3)}** outer folds over "
        f"**{n_studies}** studies and **{len(samples)}** samples."
    )

    classical: dict[str, Any] | None = None
    if CLASSICAL_PATH.is_file():
        classical = json.loads(CLASSICAL_PATH.read_text(encoding="utf-8"))
        print("loaded cached classical_baselines.json", flush=True)
    else:
        classical = run_classical(data_root=data_root, fold_pack=fold_pack, phenotypes=phenotypes)
        CLASSICAL_PATH.write_text(
            json.dumps(classical, indent=2, default=str) + "\n", encoding="utf-8"
        )

    rows.extend(classical_to_rows(classical))
    ARM_TABLE_PATH.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    plot_schema(FIG / "model_schema.png")
    plot_arm_bars(rows, FIG / "arm_bars.png")
    plot_fold_heatmap(summary, FIG / "fold_heatmap.png")
    plot_roc_from_payload(classical, FIG / "roc_tissue_ovr.png", FIG / "roc_sex.png")
    (REPORT / "analysis.md").write_text(
        write_markdown(
            rows=rows,
            summary=summary,
            classical=classical,
            n_samples=len(samples),
            n_studies_note=study_note,
        ),
        encoding="utf-8",
    )
    print(f"wrote {REPORT / 'analysis.md'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
