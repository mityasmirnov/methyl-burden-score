"""Stage 0 evaluation helpers (metrics + study-grouped splits)."""

from mbs.evaluation.annotation_slices import (
    annotation_status_counts,
    compare_hierarchical_vs_flat,
    index_annotation_summary,
    slice_metrics_from_predictions,
)
from mbs.evaluation.metrics import (
    binary_auroc_auprc,
    expected_calibration_error,
    masked_multilabel_auroc_auprc,
    metrics_by_group,
    multiclass_metrics,
    regression_metrics,
)
from mbs.evaluation.splits import (
    assert_no_study_leakage,
    build_study_grouped_split,
    grouping_key,
    partition_studies_constrained,
)

__all__ = [
    "annotation_status_counts",
    "assert_no_study_leakage",
    "binary_auroc_auprc",
    "build_study_grouped_split",
    "compare_hierarchical_vs_flat",
    "expected_calibration_error",
    "grouping_key",
    "index_annotation_summary",
    "masked_multilabel_auroc_auprc",
    "metrics_by_group",
    "multiclass_metrics",
    "partition_studies_constrained",
    "regression_metrics",
    "slice_metrics_from_predictions",
]
