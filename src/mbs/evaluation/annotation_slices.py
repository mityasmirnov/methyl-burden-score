"""Mapped vs residual evaluation slices for Milestone 6 hierarchical model."""

from __future__ import annotations

from typing import Any

import numpy as np

from mbs.batch import ANNOTATION_STATUS_NAMES
from mbs.evaluation.metrics import multiclass_metrics, regression_metrics
from mbs.training.locus_region_gene import LocusRegionGeneIndex


def annotation_status_counts(status: np.ndarray) -> dict[str, int]:
    """Count study columns per annotation status."""
    status_arr = np.asarray(status, dtype=np.int64).reshape(-1)
    out = dict.fromkeys(ANNOTATION_STATUS_NAMES, 0)
    for i, name in enumerate(ANNOTATION_STATUS_NAMES):
        out[name] = int(np.sum(status_arr == i))
    return out


def index_annotation_summary(locus_region: LocusRegionGeneIndex) -> dict[str, Any]:
    """Summarize typed vs residual topology for reports."""
    counts = annotation_status_counts(locus_region.column_annotation_status)
    return {
        "n_study_loci": locus_region.n_study_loci,
        "n_genes": locus_region.n_genes,
        "n_regions": locus_region.n_regions,
        "n_typed_edges": locus_region.n_typed_edges,
        "n_residual_cols": locus_region.n_residual_cols,
        "n_panel": locus_region.n_panel,
        "annotation_status_counts": counts,
    }


def slice_metrics_from_predictions(
    *,
    slice_name: str,
    age_true: np.ndarray | None,
    age_pred: np.ndarray | None,
    tissue_true: np.ndarray | None,
    tissue_pred: np.ndarray | None,
    sex_true: np.ndarray | None = None,
    sex_pred: np.ndarray | None = None,
    n_tissue_classes: int | None = None,
) -> dict[str, Any]:
    """Package phenotype metrics for one annotation path slice."""
    out: dict[str, Any] = {"slice": slice_name}
    if age_true is not None and age_pred is not None and len(age_true) > 0:
        out["age"] = regression_metrics(age_true, age_pred)
        out["age_n"] = len(age_true)
    if tissue_true is not None and tissue_pred is not None and len(tissue_true) > 0:
        out["tissue"] = multiclass_metrics(
            tissue_true,
            tissue_pred,
            n_classes=n_tissue_classes,
        )
        out["tissue_n"] = len(tissue_true)
    if sex_true is not None and sex_pred is not None and len(sex_true) > 0:
        out["sex"] = multiclass_metrics(sex_true, sex_pred, n_classes=2)
        out["sex_n"] = len(sex_true)
    return out


def compare_hierarchical_vs_flat(
    *,
    hierarchical_metrics: dict[str, Any],
    flat_metrics: dict[str, Any] | None,
) -> dict[str, Any]:
    """Side-by-side holdout metrics for hierarchical vs flat on the same folds."""
    comparison: dict[str, Any] = {"hierarchical": hierarchical_metrics}
    if flat_metrics is None:
        comparison["flat"] = None
        comparison["note"] = "flat metrics unavailable; compare reports manually"
        return comparison
    comparison["flat"] = flat_metrics
    deltas: dict[str, float] = {}
    for key in ("mae", "accuracy", "sex_accuracy", "loss"):
        h = hierarchical_metrics.get(key)
        f = flat_metrics.get(key)
        if isinstance(h, (int, float)) and isinstance(f, (int, float)):
            deltas[key] = float(h) - float(f)
    comparison["hierarchical_minus_flat"] = deltas
    return comparison
