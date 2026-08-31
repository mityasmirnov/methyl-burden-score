"""Classical M-value baselines for Milestone 7G (ridge / enet / HGB / PCA-SVA).

Extracted from the 7E analysis script so Hub eval can import without duplicating
fit logic. HistGradientBoosting stands in for LightGBM; PCA SVs for Bioconductor
``sva``. No metadata-only arm.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Any

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
from mbs.training.dev_cv import _phenotype_arrays
from mbs.training.features import beta_to_m_value

EPSILON = 0.001

CLASSICAL_ARMS: tuple[tuple[str, str], ...] = (
    ("C-mvalue-ridge", "ridge"),
    ("C-mvalue-enet", "enet"),
    ("C-mvalue-hgb", "hgb"),
    ("C-mvalue-sva", "sva"),
)


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


def impute_median(train: np.ndarray, test: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    imp = SimpleImputer(strategy="median")
    return imp.fit_transform(train), imp.transform(test)


def sva_residualize(
    train: np.ndarray, test: np.ndarray, n_sv: int = 10
) -> tuple[np.ndarray, np.ndarray]:
    """PCA surrogate variables (SVA-like): remove top unsupervised axes, train-only PCA."""
    n_sv = min(n_sv, train.shape[0] - 1, train.shape[1])
    pca = PCA(n_components=n_sv, svd_solver="randomized", random_state=42)
    sv_tr = pca.fit_transform(train)
    sv_te = pca.transform(test)
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


def tissue_ovr_curves(
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


def fit_eval_mvalue_fold(
    x_tr: np.ndarray,
    x_te: np.ndarray,
    ph_tr: dict[str, Any],
    ph_te: dict[str, Any],
    kind: str,
    *,
    impute: bool = True,
) -> dict[str, Any]:
    """Fit one classical kind on train; metrics on test (masked phenotypes)."""
    if impute:
        x_tr, x_te = impute_median(x_tr, x_te)
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
        age_metrics = regression_metrics(ph_te["age"][age_mt], pred)
        # Clamp exploded enet age (same policy as 7E).
        if kind == "enet" and float(age_metrics.get("mae", 0.0)) > 100.0:
            out["age"] = None
            out["age_note"] = "blanked: SGD elastic-net age MAE exploded"
        else:
            out["age"] = age_metrics
    if tissue_mt.any() and tissue_m.any():
        pred_t = tissue_model.predict(x_te[tissue_mt])
        yt = ph_te["tissue"][tissue_mt]
        tissue_valid_classes = set(ph_tr["tissue"][tissue_m].tolist())
        tissue_metrics = multiclass_metrics(yt, pred_t, valid_classes=tissue_valid_classes)
        out["tissue"] = {
            k: tissue_metrics[k]
            for k in (
                "macro_f1",
                "balanced_accuracy",
                "accuracy",
                "n_classes_scored",
                "excluded_zero_shot_test_counts",
            )
            if k in tissue_metrics
        }
        if hasattr(tissue_model, "predict_proba"):
            proba = np.asarray(tissue_model.predict_proba(x_te[tissue_mt]))
            classes = _model_classes(tissue_model)
            if classes is not None:
                out["tissue_roc"] = tissue_ovr_curves(
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


def run_classical_mvalue(
    *,
    data_root: Path,
    fold_pack: dict[str, Any],
    phenotypes: list[Any],
    max_loci: int = 65536,
    matrix_id: str = "matrix-hub-age-tissue-sex-full-v1",
) -> dict[str, Any]:
    """Run C-mvalue-* arms on frozen folds; methylation matrix only."""
    matrix_paths = matrix_store_paths(data_root / "canonical" / "matrices" / matrix_id)
    sample_index = read_sample_index(matrix_paths.sample_index_path)
    locus_index = read_locus_index(matrix_paths.locus_index_path)
    n_loci_full = int(locus_index.shape[0])
    n_cols = min(int(max_loci), n_loci_full)
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
            f"M-values on the same {n_cols}-locus prefix as the neural cascade. "
            "HistGradientBoosting stands in for LightGBM; PCA SVs for Bioconductor sva. "
            "Metadata-only omitted from ranking (7E′ leakage alarm only)."
        ),
        "folds": [],
        "roc": {},
        "arms": [a for a, _ in CLASSICAL_ARMS],
    }
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
        for arm_name, kind in CLASSICAL_ARMS:
            if kind == "sva":
                xtr_i, xte_i = impute_median(x_tr, x_te)
                xtr, xte = sva_residualize(xtr_i, xte_i)
                kind_fit = "ridge"
                extra_impute = False
            else:
                xtr, xte = x_tr, x_te
                kind_fit = kind
                extra_impute = True
            print(f"  {arm_name}…", flush=True)
            metrics = fit_eval_mvalue_fold(
                xtr, xte, ph_tr, ph_te, kind_fit, impute=extra_impute
            )
            slim = {
                k: metrics[k]
                for k in ("kind", "age", "tissue", "sex", "tissue_roc", "age_note")
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
