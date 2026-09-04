"""Shared validation checkpoint ranking (P2/P4 cascade and flat_region Stage A)."""

from __future__ import annotations

from typing import Any


def validation_rank(val_metrics: dict[str, Any]) -> tuple[float, float]:
    """Higher is better: (tissue macro-F1, -age MAE), missing -> worst."""
    f1 = val_metrics.get("macro_f1")
    if f1 is None:
        f1 = val_metrics.get("tissue_macro_f1")
    mae = val_metrics.get("mae")
    if mae is None:
        mae = val_metrics.get("age_mae")
    return (
        float(f1) if f1 is not None else -1.0,
        -float(mae) if mae is not None else -1e9,
    )


def validation_rank_age_primary(val_metrics: dict[str, Any]) -> tuple[float, float, float]:
    """Higher is better: (-age MAE, tissue macro-F1, sex AUROC). Missing -> worst."""
    mae = val_metrics.get("mae")
    if mae is None:
        mae = val_metrics.get("age_mae")
    f1 = val_metrics.get("macro_f1")
    if f1 is None:
        f1 = val_metrics.get("tissue_macro_f1")
    auroc = val_metrics.get("sex_auroc")
    if auroc is None:
        auroc = val_metrics.get("auroc")
    return (
        -float(mae) if mae is not None else -1e9,
        float(f1) if f1 is not None else -1.0,
        float(auroc) if auroc is not None else -1.0,
    )
