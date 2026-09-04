#!/usr/bin/env python3
"""Write 7G′ Stage A gene-only probe analysis.md from per_arm JSON payloads."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from mbs.annotation.manifest import write_json
from mbs.inspection.arm_glossary import render_arm_glossary_section
from mbs.inspection.comparable_metrics import (
    load_comparable_rows,
    render_comparable_ranking_section,
)
from mbs.paths import DataPaths

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REPORT = ROOT / "reports/inspection/stage0_7g_gene_only_probe"

CASCADE_ARMS = (
    "P2-G",
    "P4-G",
    "P5-G-max",
    "N-cascade-scalar-mean-max",
    "N-cascade-scalar-max-mean",
    "N-cascade-vector-mean-max",
    "N-cascade-vector-max-max",
    "N-light-gene-max",
    "N-light-gene-mean",
)
CLASSICAL_SUFFIXES = ("C-mvalue-ridge-G", "C-mvalue-enet-G", "C-mvalue-hgb-G", "C-mvalue-sva-G")
ARM_POOLING = {
    "P2-G": ("max", "max", 15),
    "P4-G": ("mean", "mean", 15),
    "P5-G-max": ("max", "max", 30),
    # Tier-1 unmatched 5-ep cells stay labeled until 16-ep promotion overwrites per_arm.
    "N-cascade-scalar-mean-max": ("mean", "max", 16),
    "N-cascade-scalar-max-mean": ("max", "mean", 16),
    "N-cascade-vector-mean-max": ("mean", "max", 16),
    "N-cascade-vector-max-max": ("max", "max", 5),
    "N-light-gene-max": ("n/a", "max", 16),
    "N-light-gene-mean": ("n/a", "mean", 16),
}
# Ceilings used only for "matched budget?" language when per_arm still has 5-ep folds.
TIER1_UNMATCHED_CEILING = 5
CASCADE_READOUTS = (
    "mbs_e2e",
    "mbs_linear_probe",
    "mbs_enet",
    "mbs_enet_nested",
    "rbs_linear_probe",
    "rbs_enet",
    "rbs_enet_nested",
)
LOCK_TRAINING_DEFAULTS = {
    "age_loss_weight": 0.3,
    "tissue_loss_weight": 3.0,
    "sex_loss_weight": 1.0,
}
CLEAR_AHEAD_DELTA = 0.03
INVALID_MBS_E2E_BANNER = (
    "> **Invalid historical `mbs_e2e`:** pre-fix runs scored train+validation+test "
    "together. Reported ~0.67–0.70 gene-only F1 is **not trustworthy**. "
    "Until rerun with `eval_split=test`, use **`mbs_linear_probe`** (~0.37) as the "
    "honest neural tissue check. **Do not lock architecture.**"
)


def _mbs_e2e_fold_valid(fold: dict[str, Any], *, arm_id: str | None = None) -> bool:
    """True when fold mbs_e2e was computed on the test split with orientation contract v2."""
    evaluations = fold.get("evaluations") or {}
    blob = evaluations.get("mbs_e2e")
    if not isinstance(blob, dict):
        return False
    if blob.get("eval_split") != "test":
        return False
    if arm_id and str(arm_id).startswith("N-light-gene"):
        manifest = fold.get("score_manifest") or {}
        if str(manifest.get("orientation_contract_version", "1")) != "2":
            return False
    return True


def _cascade_has_valid_mbs_e2e(folds: list[dict[str, Any]], *, arm_id: str | None = None) -> bool:
    if not folds:
        return False
    return all(_mbs_e2e_fold_valid(f, arm_id=arm_id) for f in folds)


def _classical_has_completed_folds(payload: dict[str, Any] | None, arm_name: str) -> bool:
    if payload is None:
        return False
    for fold in payload.get("folds") or []:
        blob = (fold.get("arms") or {}).get(arm_name) or {}
        tissue = blob.get("tissue") if isinstance(blob.get("tissue"), dict) else {}
        if tissue.get("macro_f1") is not None:
            return True
    return False


def _fmt(x: float | None) -> str:
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "—"
    return f"{x:.3f}"


def _mean_std(vals: list[float | None]) -> tuple[float | None, float | None]:
    nums = [float(v) for v in vals if v is not None and not math.isnan(v)]
    if not nums:
        return None, None
    arr = np.asarray(nums, dtype=np.float64)
    return float(arr.mean()), float(arr.std(ddof=1)) if arr.size > 1 else 0.0


def _mean_std_n(vals: list[float | None]) -> tuple[float | None, float | None, int]:
    """Mean/std plus count of non-null values (metric-specific fold count)."""
    nums = [float(v) for v in vals if v is not None and not math.isnan(v)]
    if not nums:
        return None, None, 0
    arr = np.asarray(nums, dtype=np.float64)
    std = float(arr.std(ddof=1)) if arr.size > 1 else 0.0
    return float(arr.mean()), std, int(arr.size)


def _cascade_metric_means(
    folds: list[dict[str, Any]],
    mode: str,
    *metric_keys: str,
) -> tuple[float | None, float | None]:
    path = ".".join((mode, "metrics", *metric_keys))
    per = [_metric_from_fold(f, path) for f in folds]
    return _mean_std(per)


def _cascade_metric_means_n(
    folds: list[dict[str, Any]],
    mode: str,
    *metric_keys: str,
) -> tuple[float | None, float | None, int]:
    path = ".".join((mode, "metrics", *metric_keys))
    per = [_metric_from_fold(f, path) for f in folds]
    return _mean_std_n(per)


def _classical_metric_means(
    payload: dict[str, Any],
    arm_name: str,
    *metric_keys: str,
) -> tuple[float | None, float | None]:
    per: list[float | None] = []
    for fold in payload.get("folds") or []:
        blob = (fold.get("arms") or {}).get(arm_name) or {}
        cur: Any = blob
        for key in metric_keys:
            if not isinstance(cur, dict):
                cur = None
                break
            cur = cur.get(key)
        per.append(float(cur) if cur is not None else None)
    return _mean_std(per)

def _metric_from_fold(blob: dict[str, Any], metric_path: str) -> float | None:
    parts = metric_path.split(".")
    eval_keys = (
        "mbs_e2e",
        "mbs_linear_probe",
        "mbs_enet",
        "mbs_enet_nested",
        "rbs_linear_probe",
        "rbs_enet",
        "rbs_enet_nested",
        "fusion_full",
        "fusion_mbs_direct",
    )
    if parts[0] in eval_keys:
        evaluations = blob.get("evaluations") or {}
        cur: Any = evaluations.get(parts[0])
        for key in parts[1:]:
            if not isinstance(cur, dict):
                return None
            cur = cur.get(key)
        if cur is None:
            return None
        return float(cur)
    cur: Any = blob
    for key in parts:
        if not isinstance(cur, dict):
            return None
        cur = cur.get(key)
    if cur is None:
        return None
    return float(cur)


def _cascade_tissue_f1(
    folds: list[dict[str, Any]],
    mode: str = "mbs_e2e",
) -> tuple[float | None, float | None]:
    if mode == "mbs_e2e" and not _cascade_has_valid_mbs_e2e(folds):
        return None, None
    return _cascade_metric_means(folds, mode, "tissue", "macro_f1")


def _classical_tissue_f1(payload: dict[str, Any], arm_name: str) -> tuple[float | None, float | None]:
    return _classical_metric_means(payload, arm_name, "tissue", "macro_f1")

def _gene_panel_n_cols(report_dir: Path) -> int | None:
    path = report_dir / "gene_panel_manifest.json"
    if not path.is_file():
        return None
    payload = json.loads(path.read_text(encoding="utf-8"))
    n = payload.get("n_gene_cols")
    return int(n) if n is not None else None


def load_per_arm(report_dir: Path) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    per_arm = report_dir / "per_arm"
    if not per_arm.is_dir():
        return out
    for path in sorted(per_arm.glob("*.json")):
        out[path.stem] = json.loads(path.read_text(encoding="utf-8"))
    return out


def build_lock_recommendation(
    cascade_rows: list[dict[str, Any]],
    classical_rows: list[dict[str, Any]],
    *,
    cascade_folds_by_arm: dict[str, list[dict[str, Any]]] | None = None,
    classical_payload: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Screen summary only — does **not** declare an architecture lock.

    Retained filename ``lock_recommendation.json`` for compatibility; fields
    ``locked_cascade_arm`` / ``architecture_locked`` stay null/false.
    """
    rec: dict[str, Any] = {
        "primary_metric": "mbs_e2e.metrics.tissue.macro_f1",
        "architecture_locked": False,
        "locked_cascade_arm": None,
        "best_landed_cascade_arm": None,
        "pooling_cpg": None,
        "pooling_region": None,
        "max_epochs": None,
        "cascade_clearly_ahead": None,
        "recommend_encoder_parity": False,
        "best_classical_arm": None,
        "lock_blocked_reason": (
            "ATS gene-only screen is evidence only; no cascade architecture lock. "
            "P2-G is the current reference, not a pooling lock. "
            "Next gate: matched 16-epoch promotion screen."
        ),
        "mbs_e2e_valid": False,
        "next_gate": "matched_16ep_promotion_screen",
    }
    folds_by_arm = cascade_folds_by_arm or {}
    valid_cascade_rows = [
        r
        for r in cascade_rows
        if r.get("mbs_e2e_f1") is not None
        and _cascade_has_valid_mbs_e2e(folds_by_arm.get(str(r["arm_id"]), []))
    ]
    valid_classical_rows = [
        r
        for r in classical_rows
        if r.get("tissue_f1") is not None
        and _classical_has_completed_folds(classical_payload, str(r["arm_id"]))
    ]
    if not valid_cascade_rows:
        return rec
    if not valid_classical_rows:
        rec["lock_blocked_reason"] = (
            "ATS screen incomplete (no classical -G folds); still no architecture lock. "
            "Next gate: matched 16-epoch promotion screen."
        )
        return rec

    rec["mbs_e2e_valid"] = True
    rec.update(LOCK_TRAINING_DEFAULTS)
    best_cascade = max(valid_cascade_rows, key=lambda r: r.get("mbs_e2e_f1") or -1.0)
    best_classical = max(valid_classical_rows, key=lambda r: r.get("tissue_f1") or -1.0)
    arm_id = str(best_cascade["arm_id"])
    rec["best_landed_cascade_arm"] = arm_id
    # Deliberately leave locked_cascade_arm None — do not retain a P2-G lock.
    pool = ARM_POOLING.get(arm_id)
    if pool:
        rec["pooling_cpg"], rec["pooling_region"], rec["max_epochs"] = pool
    rec["best_classical_arm"] = best_classical["arm_id"]
    c_f1 = best_cascade.get("mbs_e2e_f1")
    cl_f1 = best_classical.get("tissue_f1")
    if c_f1 is not None and cl_f1 is not None:
        rec["cascade_clearly_ahead"] = c_f1 >= cl_f1 + CLEAR_AHEAD_DELTA
        rec["recommend_encoder_parity"] = not rec["cascade_clearly_ahead"]
    return rec


def _extract_scalar(fold: dict[str, Any], *keys: str) -> float | None:
    """Drill into nested dict by sequence of keys; return float or None."""
    cur: Any = fold
    for k in keys:
        if not isinstance(cur, dict):
            return None
        cur = cur.get(k)
    try:
        return float(cur) if cur is not None else None
    except (TypeError, ValueError):
        return None


def _fold_best_epoch(fold: dict[str, Any]) -> int | None:
    """Checkpoint-selected best epoch (1-indexed when present)."""
    ckpt = fold.get("checkpoint_selection")
    if isinstance(ckpt, dict) and ckpt.get("best_epoch") is not None:
        try:
            return int(ckpt["best_epoch"])
        except (TypeError, ValueError):
            pass
    if fold.get("best_epoch") is not None:
        try:
            return int(fold["best_epoch"])
        except (TypeError, ValueError):
            pass
    diag = (
        ((fold.get("evaluations") or {}).get("mbs_e2e") or {}).get("repr_diagnostics") or {}
    )
    if diag.get("best_epoch") is not None:
        try:
            return int(diag["best_epoch"])
        except (TypeError, ValueError):
            return None
    return None


def _fold_epochs_trained(fold: dict[str, Any]) -> int | None:
    """Epochs actually run (history length / last epoch), before early stop exit."""
    history = fold.get("history")
    if isinstance(history, list) and history:
        last = history[-1]
        if isinstance(last, dict) and last.get("epoch") is not None:
            try:
                return int(last["epoch"])
            except (TypeError, ValueError):
                pass
        return len(history)
    ckpt = fold.get("checkpoint_selection")
    if isinstance(ckpt, dict):
        if ckpt.get("epochs_trained") is not None:
            try:
                return int(ckpt["epochs_trained"])
            except (TypeError, ValueError):
                pass
        vh = ckpt.get("val_history")
        if isinstance(vh, list) and vh:
            last = vh[-1]
            if isinstance(last, dict) and last.get("epoch") is not None:
                try:
                    return int(last["epoch"])
                except (TypeError, ValueError):
                    pass
            return len(vh)
    if fold.get("epochs_trained") is not None:
        try:
            return int(fold["epochs_trained"])
        except (TypeError, ValueError):
            return None
    return None


def _arm_epoch_stats(folds: list[dict[str, Any]]) -> dict[str, Any]:
    """Per-arm epoch summary for analysis.md tables."""
    bests = [e for f in folds if (e := _fold_best_epoch(f)) is not None]
    trained = [e for f in folds if (e := _fold_epochs_trained(f)) is not None]
    per_fold: list[dict[str, int | None]] = []
    for i, f in enumerate(folds):
        per_fold.append(
            {
                "fold": i,
                "best_epoch": _fold_best_epoch(f),
                "epochs_trained": _fold_epochs_trained(f),
            }
        )
    return {
        "best_epochs": bests,
        "epochs_trained": trained,
        "mean_best_epoch": float(np.mean(bests)) if bests else None,
        "mean_epochs_trained": float(np.mean(trained)) if trained else None,
        "per_fold": per_fold,
    }


def _fmt_epoch_list(vals: list[int], *, mean: float | None = None) -> str:
    if not vals:
        return "—"
    joined = ",".join(str(v) for v in vals)
    if mean is not None and len(vals) > 1:
        return f"{joined} (μ={mean:.1f})"
    return joined


def render_training_epochs_section(
    cascade_folds_by_arm: dict[str, list[dict[str, Any]]],
) -> list[str]:
    """Document epoch ceiling, epochs actually run, and best checkpoint epoch."""
    lines = [
        "## Training epochs (ceiling / ran / best)",
        "",
        "Ceiling is the configured `max_epochs` (Tier-1 screen note for N-light / "
        "mixed/vector arms). **Ran** is how many epochs the trainer completed "
        "(early stop may cut short). **Best** is the checkpoint selected by "
        "`validation_tissue_macro_f1_then_age_mae` (used for test `mbs_e2e`). "
        "Prefer actual best/ran over the ceiling label — do not stamp a hard-coded "
        "5-epoch N-light label onto longer runs.",
        "",
        "| Arm | Ceiling | Epochs ran (per fold) | Best epoch (per fold) | folds |",
        "|-----|--------:|----------------------:|----------------------:|------:|",
    ]
    any_row = False
    for arm_id in CASCADE_ARMS:
        folds = cascade_folds_by_arm.get(arm_id) or []
        if not folds:
            continue
        any_row = True
        stats = _arm_epoch_stats(folds)
        pool = ARM_POOLING.get(arm_id)
        ceiling = pool[2] if pool else "—"
        lines.append(
            f"| `{arm_id}` | {ceiling} | "
            f"{_fmt_epoch_list(stats['epochs_trained'], mean=stats['mean_epochs_trained'])} | "
            f"{_fmt_epoch_list(stats['best_epochs'], mean=stats['mean_best_epoch'])} | "
            f"{len(folds)} |"
        )
    if not any_row:
        lines.append("| — | — | — | — | 0 |")
    lines.append("")
    return lines


def _bootstrap_ci(
    vals: list[float], *, n_boot: int = 2000, alpha: float = 0.05, seed: int = 0
) -> tuple[float, float, float]:
    """Return (mean, lower, upper) bootstrap CI for a list of fold means."""
    if not vals:
        return (float("nan"),) * 3
    arr = np.asarray(vals, dtype=float)
    rng = np.random.default_rng(seed)
    boots = rng.choice(arr, size=(n_boot, len(arr)), replace=True).mean(axis=1)
    lo = float(np.percentile(boots, 100 * alpha / 2))
    hi = float(np.percentile(boots, 100 * (1 - alpha / 2)))
    return float(arr.mean()), lo, hi


def _arm_mean_ci(
    arm_payloads: dict[str, dict | None],
    arm_id: str | list[str],
    *metric_keys: str,
) -> tuple[float, float, float]:
    """Extract per-fold values then bootstrap.

    ``arm_id`` may be a single id or a list (e.g. primary + ``-s2`` seed) whose
    folds are pooled before the CI.
    """
    ids = [arm_id] if isinstance(arm_id, str) else list(arm_id)
    vals: list[float] = []
    for aid in ids:
        payload = arm_payloads.get(aid)
        if not payload:
            continue
        folds = payload.get("folds") or []
        vals.extend(v for f in folds if (v := _extract_scalar(f, *metric_keys)) is not None)
    return _bootstrap_ci(vals)


def _ablation_table_row(
    arm_id: str | list[str],
    arm_payloads: dict[str, dict | None],
    *,
    tissue_keys: tuple[str, ...],
    age_keys: tuple[str, ...],
    sex_keys: tuple[str, ...],
    linear_keys: tuple[str, ...] | None = None,
) -> str:
    t_mean, t_lo, t_hi = _arm_mean_ci(arm_payloads, arm_id, *tissue_keys)
    a_mean, a_lo, a_hi = _arm_mean_ci(arm_payloads, arm_id, *age_keys)
    s_mean, s_lo, s_hi = _arm_mean_ci(arm_payloads, arm_id, *sex_keys)

    def _fmt(m: float, lo: float, hi: float) -> str:
        if math.isnan(m):
            return "—"
        return f"{m:.3f} [{lo:.3f}–{hi:.3f}]"

    cells = [_fmt(t_mean, t_lo, t_hi), _fmt(a_mean, a_lo, a_hi), _fmt(s_mean, s_lo, s_hi)]
    if linear_keys is not None:
        l_mean, l_lo, l_hi = _arm_mean_ci(arm_payloads, arm_id, *linear_keys)
        cells.insert(1, _fmt(l_mean, l_lo, l_hi))
    return "| " + " | ".join(cells) + " |"


def _repr_diagnostics_row(arm_id: str | list[str], arm_payloads: dict[str, dict | None]) -> str:
    """Build a representation diagnostics markdown row from stored fold metrics."""
    ids = [arm_id] if isinstance(arm_id, str) else list(arm_id)
    folds: list[dict[str, Any]] = []
    for aid in ids:
        payload = arm_payloads.get(aid)
        if payload:
            folds.extend(payload.get("folds") or [])
    label = ids[0] if ids else "—"
    if not folds:
        return f"| {label} | — | — | — | — | — | — |"

    def _mean_key(*keys: str, as_int: bool = False) -> str:
        vals = [v for f in folds if (v := _extract_scalar(f, *keys)) is not None]
        vals = [v for v in vals if not (isinstance(v, float) and math.isnan(v))]
        if not vals:
            return "—"
        mean = float(np.mean(vals))
        if as_int:
            return str(int(round(mean)))
        return f"{mean:.3f}"

    score_sd = _mean_key("evaluations", "mbs_e2e", "repr_diagnostics", "gene_score_sd")
    sat_frac = _mean_key("evaluations", "mbs_e2e", "repr_diagnostics", "saturation_fraction")
    const_frac = _mean_key("evaluations", "mbs_e2e", "repr_diagnostics", "constant_score_fraction")
    corr_m = _mean_key("evaluations", "mbs_e2e", "repr_diagnostics", "corr_mean_m")
    head_l2 = _mean_key("evaluations", "mbs_e2e", "repr_diagnostics", "head_weight_l2")
    best_ep = _mean_key("evaluations", "mbs_e2e", "repr_diagnostics", "best_epoch", as_int=True)
    if best_ep == "—":
        # Fall back to fold-level checkpoint selection / best_epoch.
        be_vals = [e for f in folds if (e := _fold_best_epoch(f)) is not None]
        if be_vals:
            best_ep = str(int(round(float(np.mean(be_vals)))))
    return (
        f"| {label} | {score_sd} | {sat_frac} | {const_frac} | {corr_m} | {head_l2} | {best_ep} |"
    )


def _ablation_seed_ids(base_arm_id: str) -> list[str]:
    """Primary + optional ``-s2`` seed arm ids present under ``per_arm/``."""
    return [base_arm_id, f"{base_arm_id}-s2"]


def annotation_ablation_section(arm_payloads: dict[str, dict | None]) -> str:
    """Render annotation ablation grid (A0–A4/A7, N0–N3) with bootstrap CIs."""
    # per_arm JSONs use runner arm ids ``N-light-gene-ablation-*`` (+ ``-s2``).
    ablation_arms = [
        ("A0", "N-light-gene-ablation-m-only", "M only"),
        ("A1", "N-light-gene-ablation-m-role", "M + gene role"),
        ("A2", "N-light-gene-ablation-m-context", "M + CpG context"),
        ("A3", "N-light-gene-ablation-m-role-context", "M + role + context"),
        ("A4/A7", "N-light-gene-ablation-full", "All (regulatory zero)"),
    ]
    neg_arms = [
        ("N0", "N-light-gene-ablation-n0-obs-only", "Observed flag only"),
        ("N1", "N-light-gene-ablation-n1-anno-only", "Annotations only (no M)"),
        ("N2", "N-light-gene-ablation-n2-reg-permuted", "Reg. permuted"),
        ("N3", "N-light-gene-ablation-n3-reg-zero", "All-zero regulatory"),
    ]

    tissue_keys = ("evaluations", "mbs_e2e", "metrics", "tissue", "macro_f1")
    linear_keys = ("evaluations", "mbs_linear_probe", "metrics", "tissue", "macro_f1")
    age_keys = ("evaluations", "mbs_e2e", "metrics", "age", "mae")
    sex_keys = ("evaluations", "mbs_e2e", "metrics", "sex", "auroc")

    n_present = sum(
        1
        for _, base, _ in ablation_arms + neg_arms
        for aid in _ablation_seed_ids(base)
        if arm_payloads.get(aid)
    )

    header = [
        "## Annotation ablation grid (A0–A4, N0–N3)\n",
        "Fold 0, `mean` pooling, ≤8 epochs, **two seeds** (primary + `-s2`) pooled. "
        "Bootstrap 95% CIs over seed runs. "
        "Primary metric `mbs_e2e` tissue macro-F1; linear probe is the representation check.\n",
        f"**Payloads found:** {n_present}/18 seed runs under `per_arm/N-light-gene-ablation-*.json`.\n",
        "**Note:** A4 ≈ N2 ≈ N3 while regulatory channels are zero (cCRE/DHS/ChromHMM not "
        "on disk). **`m_only` should lead** if annotations add noise under this budget.\n",
        "| Arm | Features | Tissue e2e [95% CI] | Linear F1 [95% CI] | Age MAE (e2e) [95% CI] | Sex AUROC (e2e) [95% CI] |",
        "|-----|----------|--------------------:|-------------------:|-----------------------:|-------------------------:|",
    ]
    rows_a = [
        f"| {label} | {desc} | "
        + _ablation_table_row(
            _ablation_seed_ids(arm_id),
            arm_payloads,
            tissue_keys=tissue_keys,
            age_keys=age_keys,
            sex_keys=sex_keys,
            linear_keys=linear_keys,
        )[2:]
        for label, arm_id, desc in ablation_arms
    ]
    neg_header = [
        "",
        "### Negative controls\n",
        "| Arm | Features | Tissue e2e [95% CI] | Linear F1 [95% CI] | Age MAE (e2e) [95% CI] | Sex AUROC (e2e) [95% CI] |",
        "|-----|----------|--------------------:|-------------------:|-----------------------:|-------------------------:|",
    ]
    rows_n = [
        f"| {label} | {desc} | "
        + _ablation_table_row(
            _ablation_seed_ids(arm_id),
            arm_payloads,
            tissue_keys=tissue_keys,
            age_keys=age_keys,
            sex_keys=sex_keys,
            linear_keys=linear_keys,
        )[2:]
        for label, arm_id, desc in neg_arms
    ]

    # Short takeaway from seed-mean e2e tissue.
    ranked = []
    for label, arm_id, desc in ablation_arms + neg_arms:
        mean, _, _ = _arm_mean_ci(arm_payloads, _ablation_seed_ids(arm_id), *tissue_keys)
        if not math.isnan(mean):
            ranked.append((mean, label, desc))
    ranked.sort(reverse=True)
    takeaway = ""
    if ranked:
        best = ranked[0]
        takeaway = (
            f"\n**Takeaway:** best e2e tissue = **`{best[1]}` ({best[2]})** at "
            f"{best[0]:.3f}. "
            + (
                "`m_only` leads; gene-role/context do not help under this fold-0 budget."
                if best[1] == "A0"
                else "Unexpected leader — inspect representation / training curves."
            )
            + " Negatives `obs_only` / `anno_only` should be near chance.\n"
        )

    repr_header = [
        "",
        "### Representation diagnostics (fold 0 mean across seeds)\n",
        "Computed post-hoc from saved `scores/mbs.npy` (+ `mbs_present.npy`), "
        "checkpoint `head_state`, and Pearson r of per-sample mean MBS vs mean "
        "M-value over the gene-linked CpG panel (`sample_mean_m_gene_panel.npy`). "
        "Saturation = fraction of present scores ≤0.05 or ≥0.95; const-score = "
        "fraction of genes with SD < 1e-4 across samples.\n",
        "| Arm | Gene-score SD | Saturation frac | Const-score frac | Corr w/ mean-M | Head ‖w‖₂ | Best ep |",
        "|-----|:-------------:|:---------------:|:----------------:|:--------------:|:---------:|:-------:|",
    ]

    def _relabel(row: str, label: str) -> str:
        parts = row.split("|")
        if len(parts) < 3:
            return row
        parts[1] = f" {label} "
        return "|".join(parts)

    repr_rows_a = [
        _relabel(_repr_diagnostics_row(_ablation_seed_ids(arm_id), arm_payloads), label)
        for label, arm_id, _ in ablation_arms
    ]
    repr_rows_n = [
        _relabel(_repr_diagnostics_row(_ablation_seed_ids(arm_id), arm_payloads), label)
        for label, arm_id, _ in neg_arms
    ]

    # Brief read of diagnostics when present.
    a0_sd = _arm_mean_ci(arm_payloads, _ablation_seed_ids("N-light-gene-ablation-m-only"),
                         "evaluations", "mbs_e2e", "repr_diagnostics", "gene_score_sd")[0]
    n0_const = _arm_mean_ci(arm_payloads, _ablation_seed_ids("N-light-gene-ablation-n0-obs-only"),
                            "evaluations", "mbs_e2e", "repr_diagnostics", "constant_score_fraction")[0]
    repr_note = ""
    if not math.isnan(a0_sd) or not math.isnan(n0_const):
        bits = []
        if not math.isnan(a0_sd):
            bits.append(f"A0 gene-score SD≈{a0_sd:.3f} (non-collapsed encoder)")
        if not math.isnan(n0_const) and n0_const > 0.9:
            bits.append("N0 const-score≈1 (obs-only scores collapsed — control OK)")
        bits.append("corr(mean MBS, panel mean-M)≈0 across arms (gene MBS ≠ bulk methylation intensity)")
        bits.append("no score saturation (not stuck at 0/1)")
        repr_note = "\n**Repr read:** " + "; ".join(bits) + ".\n"

    return "\n".join(
        header
        + rows_a
        + neg_header
        + rows_n
        + [takeaway]
        + repr_header
        + repr_rows_a
        + repr_rows_n
        + [repr_note]
        + [""]
    )


def orphan_ablation_section(payload: dict[str, Any] | None) -> str:
    if payload is None:
        return "Orphan ablation arm **P2-orphan-ablation** not run.\n"
    folds = payload.get("folds") or []
    full_mean, _ = _cascade_tissue_f1(folds, "fusion_full")
    direct_mean, _ = _cascade_tissue_f1(folds, "fusion_mbs_direct")
    delta = None
    if full_mean is not None and direct_mean is not None:
        delta = full_mean - direct_mean
    lines = [
        "## Orphan RBS ablation (P2-orphan-ablation)",
        "",
        "Compare **`fusion_full`** (orphan RBS + MBS + direct) vs **`fusion_mbs_direct`** "
        "(MBS + direct only).",
        "",
        "| Mode | Mean tissue macro-F1 |",
        "|------|---------------------:|",
        f"| fusion_full | {_fmt(full_mean)} |",
        f"| fusion_mbs_direct | {_fmt(direct_mean)} |",
        f"| Δ (full − mbs_direct) | {_fmt(delta)} |",
        "",
    ]
    if delta is not None:
        if delta > 0.01:
            lines.append(
                "**Orphan RBS columns help** on the full-model fusion path (Δ > 0.01 F1)."
            )
        elif delta < -0.01:
            lines.append(
                "**Orphan RBS columns hurt** at this budget; prefer `fusion_mbs_direct` for Stage B."
            )
        else:
            lines.append(
                "**Orphan RBS effect is negligible** at this budget (|Δ| ≤ 0.01); "
                "Stage B should still report both fusion modes."
            )
    lines.append("")
    return "\n".join(lines)


def _fmt_pm(mean: float | None, std: float | None = None) -> str:
    if mean is None or (isinstance(mean, float) and math.isnan(mean)):
        return "—"
    if std is None or (isinstance(std, float) and math.isnan(std)):
        return f"{mean:.3f}"
    return f"{mean:.3f} (±{std:.3f})"


def _cascade_mode_row(folds: list[dict[str, Any]], arm_id: str, mode: str) -> dict[str, Any] | None:
    """One arm × readout row with tissue / age / sex means."""
    if mode == "mbs_e2e" and not _cascade_has_valid_mbs_e2e(folds):
        return None
    tissue, tissue_std, n_tissue = _cascade_metric_means_n(folds, mode, "tissue", "macro_f1")
    if tissue is None:
        return None
    age_mae, age_mae_std, _ = _cascade_metric_means_n(folds, mode, "age", "mae")
    age_r2, age_r2_std, _ = _cascade_metric_means_n(folds, mode, "age", "r2")
    # Nested enet age can explode under failed scaler/SGD; do not publish.
    if age_mae is not None and float(age_mae) > 100.0:
        age_mae = None
        age_mae_std = None
        age_r2 = None
        age_r2_std = None
    sex_auroc, sex_auroc_std, _ = _cascade_metric_means_n(folds, mode, "sex", "auroc")
    sex_f1, sex_f1_std, _ = _cascade_metric_means_n(folds, mode, "sex", "macro_f1")
    return {
        "arm_id": arm_id,
        "readout": mode,
        "tissue_f1": tissue,
        "tissue_f1_std": tissue_std,
        "age_mae": age_mae,
        "age_mae_std": age_mae_std,
        "age_r2": age_r2,
        "age_r2_std": age_r2_std,
        "sex_auroc": sex_auroc,
        "sex_auroc_std": sex_auroc_std,
        "sex_f1": sex_f1,
        "sex_f1_std": sex_f1_std,
        # Metric-specific: how many folds actually contribute this readout.
        "n_folds": n_tissue,
        "n_folds_total": len(folds),
    }


def _task_comparison_rows(
    cascade_folds_by_arm: dict[str, list[dict[str, Any]]],
    classical_payload: dict[str, Any] | None,
) -> list[dict[str, Any]]:
    """Side-by-side tissue / age / sex for cascade readouts and classical -G arms."""
    rows: list[dict[str, Any]] = []
    for arm_id in CASCADE_ARMS:
        folds = cascade_folds_by_arm.get(arm_id) or []
        if not folds:
            continue
        for mode in CASCADE_READOUTS:
            row = _cascade_mode_row(folds, arm_id, mode)
            if row is not None:
                rows.append(row)
    if classical_payload:
        for arm_name in CLASSICAL_SUFFIXES:
            tissue, tissue_std = _classical_tissue_f1(classical_payload, arm_name)
            if tissue is None:
                continue
            age_mae, age_mae_std = _classical_metric_means(
                classical_payload, arm_name, "age", "mae"
            )
            age_r2, age_r2_std = _classical_metric_means(classical_payload, arm_name, "age", "r2")
            sex_auroc, sex_auroc_std = _classical_metric_means(
                classical_payload, arm_name, "sex", "auroc"
            )
            sex_f1, sex_f1_std = _classical_metric_means(
                classical_payload, arm_name, "sex", "macro_f1"
            )
            n_classical = 0
            for fold in classical_payload.get("folds") or []:
                blob = (fold.get("arms") or {}).get(arm_name) or {}
                tissue_blob = blob.get("tissue") if isinstance(blob.get("tissue"), dict) else {}
                if tissue_blob.get("macro_f1") is not None:
                    n_classical += 1
            rows.append(
                {
                    "arm_id": arm_name,
                    "readout": "classical",
                    "tissue_f1": tissue,
                    "tissue_f1_std": tissue_std,
                    "age_mae": age_mae,
                    "age_mae_std": age_mae_std,
                    "age_r2": age_r2,
                    "age_r2_std": age_r2_std,
                    "sex_auroc": sex_auroc,
                    "sex_auroc_std": sex_auroc_std,
                    "sex_f1": sex_f1,
                    "sex_f1_std": sex_f1_std,
                    "n_folds": n_classical,
                    "n_folds_total": len(classical_payload.get("folds") or []),
                }
            )
    return rows


def render_task_comparison_section(rows: list[dict[str, Any]]) -> list[str]:
    """Single table: tissue F1, age MAE/R², sex AUROC/F1."""
    lines = [
        "## Task comparison (tissue / age / sex)",
        "",
        "Same **`explicit_only`** gene-linked panel and outer **test** folds. "
        "Compare rows as alternative **readouts** of one encoder (`mbs_e2e` / "
        "`mbs_linear_probe` / `mbs_enet` / `rbs_*`) versus classical models on raw "
        "CpG M-values. Prefer `mbs_e2e` sex AUROC when present (proba path); "
        "`rbs_*` diagnose loss before vs after gene pooling. "
        "Classical enet age uses Huber SGD elastic-net (year-scale target + "
        "eta0=1e-4); unscaled squared-error SGD exploded on this panel. "
        "Horvath-style clocks are not in this table.",
        "",
        "| Arm | Readout | Tissue F1 | Age MAE | Age R² | Sex AUROC | Sex F1 | folds |",
        "|-----|---------|----------:|--------:|-------:|----------:|-------:|------:|",
    ]
    ranked = sorted(
        rows,
        key=lambda r: r.get("tissue_f1") if r.get("tissue_f1") is not None else -1.0,
        reverse=True,
    )
    for row in ranked:
        lines.append(
            f"| `{row['arm_id']}` | `{row['readout']}` | "
            f"{_fmt_pm(row.get('tissue_f1'), row.get('tissue_f1_std'))} | "
            f"{_fmt_pm(row.get('age_mae'), row.get('age_mae_std'))} | "
            f"{_fmt_pm(row.get('age_r2'), row.get('age_r2_std'))} | "
            f"{_fmt_pm(row.get('sex_auroc'), row.get('sex_auroc_std'))} | "
            f"{_fmt_pm(row.get('sex_f1'), row.get('sex_f1_std'))} | "
            f"{row.get('n_folds', 0)} |"
        )
    lines.extend(
        [
            "",
            "**Readouts:** `mbs_e2e` = jointly trained neural heads on MBS (ATS screen primary); "
            "`mbs_linear_probe` / `mbs_enet` = new sklearn heads on the **same frozen MBS**; "
            "`rbs_linear_probe` / `rbs_enet` = frozen **gene-linked RBS** (pre–gene-pool); "
            "`classical` = sklearn on gene-linked CpG M-values (no encoder). "
            "**folds** = number of folds that actually contain that readout "
            "(±0.000 with folds=1 means a single fold, not three identical scores).",
            "",
        ]
    )
    return lines


def _pareto_front(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Non-dominated on tissue F1 (↑), age MAE (↓), sex AUROC (↑). Prefer mbs_e2e."""
    candidates = [
        r
        for r in rows
        if r.get("readout") in ("mbs_e2e", "classical")
        and r.get("tissue_f1") is not None
        and r.get("age_mae") is not None
    ]
    front: list[dict[str, Any]] = []
    for a in candidates:
        dominated = False
        for b in candidates:
            if a is b:
                continue
            better_or_eq = (
                (b["tissue_f1"] >= a["tissue_f1"])
                and (b["age_mae"] <= a["age_mae"])
                and (
                    (b.get("sex_auroc") or 0.0) >= (a.get("sex_auroc") or 0.0)
                )
            )
            strictly = (
                b["tissue_f1"] > a["tissue_f1"]
                or b["age_mae"] < a["age_mae"]
                or (b.get("sex_auroc") or 0.0) > (a.get("sex_auroc") or 0.0)
            )
            if better_or_eq and strictly:
                dominated = True
                break
        if not dominated:
            front.append(a)
    return sorted(front, key=lambda r: float(r["tissue_f1"]), reverse=True)


def render_pareto_section(rows: list[dict[str, Any]]) -> list[str]:
    front = _pareto_front(rows)
    lines = [
        "## Three-task Pareto (`mbs_e2e` + classical)",
        "",
        "Non-dominated on tissue macro-F1 (↑), age MAE (↓), sex AUROC (↑). "
        "Do **not** pick a winner on tissue alone.",
        "",
        "| Arm | Readout | Tissue F1 | Age MAE | Sex AUROC |",
        "|-----|---------|----------:|--------:|----------:|",
    ]
    for row in front:
        lines.append(
            f"| `{row['arm_id']}` | `{row['readout']}` | "
            f"{_fmt(row.get('tissue_f1'))} | {_fmt(row.get('age_mae'))} | "
            f"{_fmt(row.get('sex_auroc'))} |"
        )
    if not front:
        lines.append("| — | — | — | — | — |")
    lines.extend(["", ""])
    return lines


def render_architecture_qa(
    rows: list[dict[str, Any]],
    *,
    cascade_folds_by_arm: dict[str, list[dict[str, Any]]] | None = None,
) -> list[str]:
    """Answer the seven Stage A architecture questions from available rows."""
    folds_by_arm = cascade_folds_by_arm or {}

    def _e2e(arm: str) -> dict[str, Any] | None:
        for r in rows:
            if r.get("arm_id") == arm and r.get("readout") == "mbs_e2e":
                return r
        return None

    def _mode(arm: str, mode: str) -> dict[str, Any] | None:
        for r in rows:
            if r.get("arm_id") == arm and r.get("readout") == mode:
                return r
        return None

    p2 = _e2e("P2-G")
    p4 = _e2e("P4-G")
    mean_max = _e2e("N-cascade-scalar-mean-max")
    max_mean = _e2e("N-cascade-scalar-max-mean")
    vec = _e2e("N-cascade-vector-mean-max") or _e2e("N-cascade-vector-max-max")
    light = _e2e("N-light-gene-max") or _e2e("N-light-gene-mean")
    # Prefer locked P2 if present; else the screen arm with the clearest
    # RBS-vs-MBS age gap (gene-pooling diagnostic), falling back to any RBS row.
    rbs_arm = None
    rbs = None
    mbs_probe = None
    candidates: list[tuple[str, dict[str, Any], dict[str, Any], float]] = []
    for arm in (
        "P2-G",
        "N-cascade-scalar-max-mean",
        "N-cascade-scalar-mean-max",
        "N-cascade-vector-mean-max",
        "N-cascade-vector-max-max",
        "P4-G",
    ):
        # Prefer linear for vector RBS (fixed enet is an unsuitable sparse readout).
        if "vector" in arm:
            cand = _mode(arm, "rbs_linear_probe") or _mode(arm, "rbs_enet")
        else:
            cand = _mode(arm, "rbs_enet") or _mode(arm, "rbs_linear_probe")
        if cand is None:
            continue
        probe = _mode(arm, "mbs_enet") or _mode(arm, "mbs_linear_probe")
        if probe is None:
            continue
        rbs_age = float(cand.get("age_mae") or 99.0)
        mbs_age = float(probe.get("age_mae") or 99.0)
        gap = mbs_age - rbs_age  # >0 ⇒ gene pool worsens age
        candidates.append((arm, cand, probe, gap))
    if candidates:
        # Prefer P2 if present; else largest positive age gap.
        p2_hit = next((c for c in candidates if c[0] == "P2-G"), None)
        chosen = p2_hit or max(candidates, key=lambda c: c[3])
        rbs_arm, rbs, mbs_probe, _ = chosen
    classical = next((r for r in rows if r.get("arm_id") == "C-mvalue-enet-G"), None)

    def _arm_ceiling(arm_id: str | None) -> int | None:
        if arm_id is None:
            return None
        pool = ARM_POOLING.get(arm_id)
        return int(pool[2]) if pool else None

    def _fold_ceilings(arm_id: str) -> list[int]:
        """Actual configured ceilings from fold checkpoint_selection when present."""
        out: list[int] = []
        for fold in folds_by_arm.get(arm_id) or []:
            ckpt = fold.get("checkpoint_selection") or {}
            if isinstance(ckpt, dict) and ckpt.get("max_epochs") is not None:
                try:
                    out.append(int(ckpt["max_epochs"]))
                    continue
                except (TypeError, ValueError):
                    pass
            trained = _fold_epochs_trained(fold)
            if trained is not None:
                out.append(int(trained))
        return out

    def _budgets_matched(arm_a: str | None, arm_b: str | None) -> bool:
        """True when both arms ran near the same epoch ceiling (15≈16 OK)."""
        if arm_a is None or arm_b is None:
            return False
        ca = _fold_ceilings(arm_a)
        cb = _fold_ceilings(arm_b)
        if not ca:
            ca_fallback = _arm_ceiling(arm_a)
            ca = [ca_fallback] if ca_fallback is not None else []
        if not cb:
            cb_fallback = _arm_ceiling(arm_b)
            cb = [cb_fallback] if cb_fallback is not None else []
        if not ca or not cb:
            return False
        # Any fold still at the Tier-1 5-ep ceiling is unmatched vs P2/P4.
        if any(c <= TIER1_UNMATCHED_CEILING for c in ca + cb):
            return False
        return abs(max(ca) - max(cb)) <= 1

    def _cmp_pool(level: str) -> str:
        if level == "cpg":
            a, b = mean_max, p2  # mean vs max at cpg (region fixed max)
            label_a, label_b = "mean-max", "max-max"
            arm_a, arm_b = "N-cascade-scalar-mean-max", "P2-G"
            if mean_max is None or p2 is None:
                a, b = p4, max_mean
                label_a, label_b = "mean-mean", "max-mean"
                arm_a, arm_b = "P4-G", "N-cascade-scalar-max-mean"
        else:
            a, b = max_mean, p2  # mean vs max at gene (cpg fixed max)
            label_a, label_b = "max-mean", "max-max"
            arm_a, arm_b = "N-cascade-scalar-max-mean", "P2-G"
            if max_mean is None or p2 is None:
                a, b = p4, mean_max
                label_a, label_b = "mean-mean", "mean-max"
                arm_a, arm_b = "P4-G", "N-cascade-scalar-mean-max"
        if a is None or b is None:
            return f"Insufficient arms for {level}-level comparison yet."
        nums = (
            f"`{label_a}` tissue F1={_fmt(a.get('tissue_f1'))} vs "
            f"`{label_b}` {_fmt(b.get('tissue_f1'))}; "
            f"age MAE {_fmt(a.get('age_mae'))} vs {_fmt(b.get('age_mae'))}."
        )
        if not _budgets_matched(arm_a, arm_b):
            return (
                f"{nums} **Unmatched epoch budgets — no pooling lock.** "
                f"`P2-G` is the current reference (best 15-ep scalar result), not a "
                f"resolved pooling winner; promote mixed cells to a matched ceiling first."
            )
        better = a if (a.get("tissue_f1") or 0) >= (b.get("tissue_f1") or 0) else b
        return (
            f"{nums} Prefer **`{better['arm_id']}`** on this matched slice (check Pareto)."
        )

    q3 = "Pending RBS diagnostic."
    if rbs is not None and mbs_probe is not None:
        q3 = (
            f"`{rbs_arm}` `rbs_*` tissue F1={_fmt(rbs.get('tissue_f1'))}, "
            f"age MAE={_fmt(rbs.get('age_mae'))}, sex AUROC={_fmt(rbs.get('sex_auroc'))}; "
            f"same-arm MBS probe tissue={_fmt(mbs_probe.get('tissue_f1'))}, "
            f"age={_fmt(mbs_probe.get('age_mae'))}, sex={_fmt(mbs_probe.get('sex_auroc'))}. "
            + (
                "Gene pooling is near-neutral on tissue; **age/sex often better on RBS** "
                "(pre–gene-pool), so some phenotype signal is lost at region→gene."
                if (rbs.get("age_mae") or 99) + 0.5 < (mbs_probe.get("age_mae") or 0)
                or (rbs.get("sex_auroc") or 0) > (mbs_probe.get("sex_auroc") or 0) + 0.02
                else "Loss is not clearly at gene pooling on this arm; check scalar RBS vs raw CpG."
            )
        )
    elif rbs is not None:
        q3 = (
            f"`{rbs_arm}` `rbs_*` tissue F1={_fmt(rbs.get('tissue_f1'))}, "
            f"age MAE={_fmt(rbs.get('age_mae'))} (no paired MBS probe row)."
        )
    if classical is not None and rbs is not None:
        q3 += (
            f" Classical enet age MAE={_fmt(classical.get('age_mae'))} remains the age ceiling."
        )

    q5 = "Pending one-hop arms."
    if light is not None and p2 is not None:
        q5 = (
            f"One-hop `{light['arm_id']}` tissue={_fmt(light.get('tissue_f1'))} / "
            f"age={_fmt(light.get('age_mae'))} vs P2-G "
            f"{_fmt(p2.get('tissue_f1'))} / {_fmt(p2.get('age_mae'))}."
        )

    q6 = (
        "Gene aggregation still trails classical on age/sex; one scalar MBS/gene is "
        "**not yet adequate** unless a screen arm closes the gap."
        if classical is not None
        else "Pending classical comparison."
    )

    lines = [
        "## Architecture questions (Stage A screen)",
        "",
        f"1. **CpG → region pool (mean vs max):** {_cmp_pool('cpg')}",
        f"2. **Region → gene pool (mean vs max):** {_cmp_pool('region')}",
        f"3. **Does scalar RBS discard information?** Vector arm "
        f"`{_fmt(vec.get('tissue_f1') if vec else None)}` tissue vs P2 "
        f"`{_fmt(p2.get('tissue_f1') if p2 else None)}`; "
        + (
            "if vector does not beat scalar on age/sex, bottleneck is elsewhere."
            if vec is not None
            else "vector arms pending."
        ),
        f"4. **Gene pooling vs RBS:** {q3}",
        f"5. **One-hop vs cascade:** {q5}",
        f"6. **One-scalar-per-gene bottleneck:** {q6}",
        "7. **Best performance/compute:** Prefer landed P2/P4 (15 ep) as the current "
        "ATS **reference**, not a pooling lock. Do **not** promote unmatched Tier-1 "
        "(5 ep) arms against 15-ep P2. Next gate is the **matched 16-epoch promotion "
        "screen**; age-primary seed-mask waits on those decision rules.",
        "",
    ]
    return lines


def render_rbs_frozen_section(task_rows: list[dict[str, Any]]) -> list[str]:
    """Post-hoc / inline RBS readouts for screen cascade arms."""
    screen_arms = (
        "N-cascade-scalar-mean-max",
        "N-cascade-scalar-max-mean",
        "N-cascade-vector-mean-max",
        "N-cascade-vector-max-max",
    )
    by_arm: dict[str, dict[str, dict[str, Any]]] = {a: {} for a in screen_arms}
    for row in task_rows:
        arm = str(row.get("arm_id") or "")
        mode = str(row.get("readout") or "")
        if arm in by_arm and mode in ("rbs_enet", "rbs_linear_probe"):
            by_arm[arm][mode] = row
    if not any(by_arm[a] for a in screen_arms):
        return []
    lines = [
        "",
        "### RBS frozen readouts (screen cascade — `rbs_enet` / `rbs_linear`)",
        "",
        "| Arm | `rbs_enet` tissue F1 | `rbs_enet` age MAE | `rbs_enet` sex AUROC | "
        "`rbs_linear` tissue F1 | `rbs_linear` age MAE | folds |",
        "|-----|---------------------:|-------------------:|---------------------:"
        "|-----------------------:|---------------------:|------:|",
    ]
    for arm in screen_arms:
        enet = by_arm[arm].get("rbs_enet")
        lin = by_arm[arm].get("rbs_linear_probe")
        if enet is None and lin is None:
            continue
        n = (enet or lin or {}).get("n_folds") or 0
        lines.append(
            f"| `{arm}` | "
            f"{_fmt_pm((enet or {}).get('tissue_f1'), (enet or {}).get('tissue_f1_std'))} | "
            f"{_fmt_pm((enet or {}).get('age_mae'), (enet or {}).get('age_mae_std'))} | "
            f"{_fmt_pm((enet or {}).get('sex_auroc'), (enet or {}).get('sex_auroc_std'))} | "
            f"{_fmt_pm((lin or {}).get('tissue_f1'), (lin or {}).get('tissue_f1_std'))} | "
            f"{_fmt_pm((lin or {}).get('age_mae'), (lin or {}).get('age_mae_std'))} | "
            f"{n} |"
        )
    lines.extend(
        [
            "",
            "`rbs_enet` via `scripts/eval_mbs_enet_from_scores.py --which rbs` on saved "
            "`all_gene_rbs.zarr` (13 212 regions; no encoder retrain). Fixed "
            "`alpha=0.1` / `l1_ratio=0.5` **without** train-fold standardization is "
            "**diagnostic only**. Scalar arms: enet ≈/≥ linear on tissue and improves "
            "age. Vector arms: age collapses under that fixed enet while sex stays "
            "nearly unchanged — that is an over-strong / unscaled sparse penalty, "
            "**not** evidence the vector RBS representation is weak. Prefer "
            "`rbs_linear_probe` (and nested `rbs_enet_nested` once available) for "
            "vector RBS. P2-G `rbs_enet` not run (folds 1–2 lack `all_gene_rbs.zarr`).",
            "",
        ]
    )
    return lines


def _interpretation_section() -> list[str]:
    """Durable scientific interpretation (must survive report regeneration)."""
    return [
        "## Interpretation",
        "",
        "### N-light is not collapsed",
        "",
        "N-light mean improves from age MAE ~23.08 end-to-end to ~11.55 with a refitted "
        "linear head. The encoder contains information; the native optimisation / readout "
        "is poor.",
        "",
        "### Annotation graph is populated; network input is weak",
        "",
        "Audit: 51,375 unique CpGs; 57,430 locus–gene edges; 2,646 genes; 5,718 "
        "multi-gene CpGs; zero `other_gene` edges; all CpG-context categories populated. "
        "The short raw-concatenation ablation does **not** show annotations are "
        "uninformative — it shows **raw concatenation** hurts a short, tissue-primary, "
        "mean-pooling run. Implementation weaknesses (`gather_flat_region_features` / "
        "`FlatDeepSetRegion`):",
        "",
        "- M-value + six gene-role one-hots + seven context one-hots + six regulatory "
        "slots + flags → one 24-d input;",
        "- all six regulatory channels currently zero;",
        "- `observed` effectively always one (unobserved edges dropped before encode);",
        "- `gene_role_present` effectively constant on the Stage A graph;",
        "- raw unnormalized M mixed with 0/1 annotations;",
        "- global max/mean pool lets promoter and body cancel;",
        "- fold-fitted robust-z fitted in `loop.py` but unused by the flat-region path.",
        "",
        "Preferred one-hop (document only; do not train in this gate): gated embeddings "
        "`h_i = φ_M(z_i) + α_R E_R(role) + α_C E_C(context)` with `α` near 0, then "
        "**role-stratified** pools (promoter-core / proximal / 5′ / body / 3′) → gene "
        "embedding → scalar MBS.",
        "",
        "### Vector vs scalar",
        "",
        "Fair five-epoch mean→max: scalar tissue F1 ~0.331 vs vector ~0.337. Do **not** "
        "compare five-epoch vector to fifteen-epoch P2. Vector RBS **linear** probe "
        "before gene pool: age MAE ~10.46, sex AUROC ~0.842 — CpG→region works; failure "
        "is later (elementwise max/mean, no typed output channel, one scalar MBS). "
        "Fixed **`rbs_enet` on vector mean→max collapses age** (MAE ~19.67) while sex "
        "holds (~0.837 vs linear 0.842). That pattern is characteristic of an overly "
        "strong or poorly calibrated sparse penalty on a distributed age signal — "
        "**not** evidence that the vector RBS representation is weak. Prefer "
        "`rbs_linear_probe` / nested enet for vector diagnostics. Raising LR is not "
        "the first response.",
        "",
        "### RBS → MBS and typed pooling",
        "",
        "Cascade already adds a region-type embedding before RBS; the scalar path then "
        "pools only scalars (type discarded) and the vector path can still mix roles. "
        "CPU ablation **R0–R5** "
        "(`reports/inspection/stage0_7g_gene_only_probe/typed_rbs_pooling/`) finds typed "
        "max/mean (R1–R3) improve age MAE by ~3–4 y vs untyped R0 on vector-mean-max "
        "RBS, but the within-gene **role shuffle control does not collapse** — so the "
        "gain is not yet proof of biological role identity (extra channels / capacity "
        "may explain it). This diagnostic does **not** decide Stage B. A neural typed "
        "aggregator remains a pooling follow-up only if typed arms beat R0 on age "
        "**and** the shuffle control collapses.",
        "",
        "### Next real gate",
        "",
        "**Matched 16-epoch promotion screen** (one-hop max/mean, scalar mixed pools, "
        "vector mean→max) before declaring pooling winners or starting age-primary "
        "seed-mask. `P2-G` is the current reference, **not** a pooling lock. "
        "Fold-selected-panel Stage B stays blocked.",
        "",
    ]


def write_analysis(report_dir: Path, *, lock: dict[str, Any], paths: DataPaths | None = None) -> None:
    arms = load_per_arm(report_dir)
    classical_payload = arms.get("C-mvalue-classical-G")
    classical_rows: list[dict[str, Any]] = []
    cascade_rows: list[dict[str, Any]] = []
    cascade_folds_by_arm: dict[str, list[dict[str, Any]]] = {}
    if classical_payload:
        for arm_name in CLASSICAL_SUFFIXES:
            f1, f1_std = _classical_tissue_f1(classical_payload, arm_name)
            age_mae, _ = _classical_metric_means(classical_payload, arm_name, "age", "mae")
            age_r2, _ = _classical_metric_means(classical_payload, arm_name, "age", "r2")
            sex_auroc, _ = _classical_metric_means(classical_payload, arm_name, "sex", "auroc")
            classical_rows.append(
                {
                    "arm_id": arm_name,
                    "tissue_f1": f1,
                    "tissue_f1_std": f1_std,
                    "age_mae": age_mae,
                    "age_r2": age_r2,
                    "sex_auroc": sex_auroc,
                }
            )

    for arm_id in CASCADE_ARMS:
        payload = arms.get(arm_id)
        if payload is None:
            continue
        folds = payload.get("folds") or []
        cascade_folds_by_arm[arm_id] = folds
        e2e, e2e_std, n_e2e = _cascade_metric_means_n(folds, "mbs_e2e", "tissue", "macro_f1")
        if not _cascade_has_valid_mbs_e2e(folds, arm_id=arm_id):
            e2e, e2e_std, n_e2e = None, None, 0
        probe, _, n_probe = _cascade_metric_means_n(folds, "mbs_linear_probe", "tissue", "macro_f1")
        enet, enet_std, n_enet = _cascade_metric_means_n(folds, "mbs_enet", "tissue", "macro_f1")
        age_mae, _, _ = _cascade_metric_means_n(folds, "mbs_e2e", "age", "mae")
        age_r2, _, _ = _cascade_metric_means_n(folds, "mbs_e2e", "age", "r2")
        sex_auroc, _, _ = _cascade_metric_means_n(folds, "mbs_linear_probe", "sex", "auroc")
        if sex_auroc is None:
            sex_auroc, _, _ = _cascade_metric_means_n(folds, "mbs_enet", "sex", "auroc")
        contaminated_e2e = None
        if e2e is None and folds:
            contaminated_e2e, _ = _cascade_metric_means(folds, "mbs_e2e", "tissue", "macro_f1")
        cascade_rows.append(
            {
                "arm_id": arm_id,
                "mbs_e2e_f1": e2e,
                "mbs_e2e_std": e2e_std,
                "mbs_e2e_contaminated_f1": contaminated_e2e,
                "mbs_e2e_valid": _cascade_has_valid_mbs_e2e(folds, arm_id=arm_id),
                "mbs_linear_probe_f1": probe,
                "mbs_enet_f1": enet,
                "mbs_enet_std": enet_std,
                "age_mae": age_mae if e2e is not None else None,
                "age_r2": age_r2 if e2e is not None else None,
                "sex_auroc": sex_auroc if e2e is not None else None,
                "n_folds": len(folds),
                "n_folds_e2e": n_e2e,
                "n_folds_probe": n_probe,
                "n_folds_enet": n_enet,
            }
        )

    if not lock:
        lock = build_lock_recommendation(
            cascade_rows,
            classical_rows,
            cascade_folds_by_arm=cascade_folds_by_arm,
            classical_payload=classical_payload,
        )
    write_json(report_dir / "lock_recommendation.json", lock)

    glossary_ids: list[str] = list(CASCADE_ARMS) + list(CLASSICAL_SUFFIXES) + [
        "C-mvalue-classical-G",
        "P2-orphan-ablation",
        "N-cascade-scalar-max-max",
        "N-cascade-scalar-mean-mean",
        "rbs_linear_probe",
        "rbs_enet",
        "C-mvalue-enetS",
        "N-cascade-S",
        "N-light-type",
        "N-mbs-posthoc-full-fusion",
        "N-mbs-posthoc-mbs-direct",
    ]
    if lock.get("best_landed_cascade_arm"):
        glossary_ids.insert(0, str(lock["best_landed_cascade_arm"]))

    lines = [
        "# 7G′ Stage A — gene-only MBS architecture selection",
        "",
        "Primary metric: **`mbs_e2e`** tissue macro-F1 (end-to-end MBS heads; not late fusion).",
        "Classical comparator: **`C-mvalue-*-G`** on identical `gene_cols`.",
        "",
    ]
    if any(not r.get("mbs_e2e_valid") for r in cascade_rows if r.get("n_folds")):
        lines.extend([INVALID_MBS_E2E_BANNER, ""])
    lines.extend(
        render_arm_glossary_section(
            glossary_ids,
            extra_eval_modes=(
                "mbs_e2e",
                "mbs_linear_probe",
                "mbs_enet",
                "rbs_linear_probe",
                "rbs_enet",
                "fusion_full",
                "fusion_mbs_direct",
            ),
        )
    )
    if paths is not None:
        classical_path = (
            report_dir.parent / "stage0_7g_methylation_eval" / "classical_baselines.json"
        )
        comp_rows = load_comparable_rows(
            paths.artifact_root,
            classical_baselines_path=classical_path,
            gene_classical_path=report_dir / "per_arm" / "C-mvalue-classical-G.json",
        )
        write_json(report_dir / "comparable_ranking.json", {"rows": comp_rows})
        lines.extend(render_comparable_ranking_section(comp_rows))
    task_rows = _task_comparison_rows(cascade_folds_by_arm, classical_payload)
    write_json(report_dir / "task_comparison.json", {"rows": task_rows})
    lines.extend(render_task_comparison_section(task_rows))
    lines.extend(render_pareto_section(task_rows))
    lines.extend(render_architecture_qa(task_rows, cascade_folds_by_arm=cascade_folds_by_arm))
    lines.extend(render_training_epochs_section(cascade_folds_by_arm))
    lines.extend(
        [
        "## Cascade arms (gene-linked CpGs only)",
        "",
        "Primary **`mbs_e2e`** (test split only); **`mbs_linear_probe`** and **`mbs_enet`** "
        "are readouts of the **same frozen MBS**; **`rbs_linear_probe`** / **`rbs_enet`** "
        "use gene-linked RBS (`all_gene_rbs.zarr`). "
        "Contaminated pre-fix **`mbs_e2e`** shown as *invalid*. "
        "**Best ep** = checkpoint epoch used for test eval; **ran** = epochs completed.",
        "",
        "| Arm | mbs_e2e F1 | linear probe F1 | mbs_enet F1 | age MAE (e2e) | sex AUROC (probe) | best ep | ran | folds |",
        "|-----|-----------:|----------------:|------------:|--------------:|------------------:|--------:|----:|------:|",
        ]
    )
    for row in sorted(
        cascade_rows,
        key=lambda r: r.get("mbs_e2e_f1") or r.get("mbs_linear_probe_f1") or -1.0,
        reverse=True,
    ):
        e2e_disp = _fmt(row.get("mbs_e2e_f1"))
        if not row.get("mbs_e2e_valid") and str(row.get("arm_id", "")).startswith("N-light-gene"):
            e2e_disp = f"*invalid (pre-fix orient)* {e2e_disp}"
        elif row.get("mbs_e2e_std") is not None and row.get("mbs_e2e_f1") is not None:
            e2e_disp = f"{e2e_disp} (±{_fmt(row.get('mbs_e2e_std'))})"
        elif row.get("mbs_e2e_f1") is None:
            e2e_disp = "—"
        enet_disp = _fmt(row.get("mbs_enet_f1"))
        if row.get("mbs_enet_std") is not None and row.get("mbs_enet_f1") is not None:
            enet_disp = f"{enet_disp} (±{_fmt(row.get('mbs_enet_std'))})"
        n_enet = row.get("n_folds_enet")
        if n_enet is not None and row.get("mbs_enet_f1") is not None:
            n_total = row.get("n_folds") or 0
            if int(n_enet) < int(n_total):
                enet_disp = f"{enet_disp} [{n_enet}/{n_total}]"
        folds = cascade_folds_by_arm.get(str(row["arm_id"])) or []
        ep = _arm_epoch_stats(folds)
        best_disp = _fmt_epoch_list(ep["best_epochs"], mean=ep["mean_best_epoch"])
        ran_disp = _fmt_epoch_list(ep["epochs_trained"], mean=ep["mean_epochs_trained"])
        lines.append(
            f"| {row['arm_id']} | {e2e_disp} | "
            f"{_fmt(row.get('mbs_linear_probe_f1'))} | {enet_disp} | "
            f"{_fmt(row.get('age_mae'))} | "
            f"{_fmt(row.get('sex_auroc'))} | {best_disp} | {ran_disp} | "
            f"{row.get('n_folds', 0)} |"
        )
    lines.extend(render_rbs_frozen_section(task_rows))
    n_gene_cols = _gene_panel_n_cols(report_dir)
    panel_label = f"**{n_gene_cols:,} gene-linked CpGs**" if n_gene_cols else "gene-linked CpGs"
    lines.extend(
        [
            "",
            "## Classical arms (-G panel)",
            "",
            f"Same {panel_label} as neural arms (`explicit_only`): ridge, elastic-net, HGB, PCA-SVA+ridge.",
            "",
            "| Arm | tissue F1 | age MAE | age R² | sex AUROC |",
            "|-----|----------:|--------:|-------:|----------:|",
        ]
    )
    for row in sorted(classical_rows, key=lambda r: r.get("tissue_f1") or -1.0, reverse=True):
        lines.append(
            f"| {row['arm_id']} | {_fmt(row.get('tissue_f1'))} "
            f"(±{_fmt(row.get('tissue_f1_std'))}) | "
            f"{_fmt(row.get('age_mae'))} | {_fmt(row.get('age_r2'))} | "
            f"{_fmt(row.get('sex_auroc'))} |"
        )

    lines.extend(
        [
            "",
            "## Screen status (no architecture lock)",
            "",
            "The ATS gene-only screen is **evidence**, not an architecture decision. "
            "**No cascade topology is locked.** Fold-selected-panel Stage B is **not** "
            "the next gate.",
            "",
        ]
    )
    best = lock.get("best_landed_cascade_arm")
    if best:
        lines.extend(
            [
                f"- **Best landed ATS cascade row:** `{best}` "
                f"(pooling `{lock.get('pooling_cpg')}` / `{lock.get('pooling_region')}`; "
                f"configured ceiling {lock.get('max_epochs')} ep — see Training epochs for actual best/ran)",
                f"- **Best classical (-G):** `{lock.get('best_classical_arm')}`",
                f"- **Cascade clearly ahead (≥{CLEAR_AHEAD_DELTA} tissue F1):** {lock.get('cascade_clearly_ahead')}",
                f"- **Architecture locked:** `{lock.get('architecture_locked', False)}`",
                "",
                "Caveats: not ≥0.03 ahead of classical; tissue-primary loss; "
                "`P2-G` is the **current reference**, not a pooling lock. Unmatched "
                "Tier-1 (5 ep) cells must not be compared to 15-ep P2 as if budgets "
                "matched — run the 16-epoch promotion screen first.",
                "",
            ]
        )
    else:
        reason = lock.get("lock_blocked_reason") or "insufficient evidence"
        lines.extend(
            [
                f"- **Status:** {reason}",
                "",
            ]
        )

    lines.append(orphan_ablation_section(arms.get("P2-orphan-ablation")))
    lines.append(annotation_ablation_section(arms))
    lines.extend(_interpretation_section())
    lines.extend(
        [
            "## Parallel / follow-on work",
            "",
            "- **ATS Stage A Tier-1 screen + annotation ablations:** complete. "
            "Freeze **`P2-G` as current reference, not a pooling lock.**",
            "- **Matched 16-epoch promotion screen (current GPU gate):** "
            "**`N-light-gene-max`** tissue e2e **0.336** (below P2 — do not rerun); "
            "**`N-light-gene-mean`** tissue e2e **0.378** (**within ~0.03 of P2** — "
            "`one_hop_mean_near_p2` fired); **`N-cascade-scalar-mean-max` 16-ep** "
            "done 3/3 (mean e2e ≈0.346). **Now:** scalar **max→mean** 16-ep "
            "(fold 0 done; fold 1+ training), then vector mean→max ×3 — "
            "[`milestone-7g-prime-16ep-promotion.md`]"
            "(../../../docs/plans/milestone-7g-prime-16ep-promotion.md).",
            "- **Post-hoc CPU enet:** light max/mean 16-ep already have nested "
            "`mbs_enet`. Launching / filling fixed+nested **mbs/rbs** enet on "
            "`scalar-mean-max-16ep` and available `scalar-max-mean-16ep` folds "
            "(`scratch/logs/16ep_posthoc_enet.log`; unit `mbs-16ep-enet`). "
            "Fixed `rbs_enet` remains diagnostic-only.",
            "- **CPU typed-RBS ablation (R0–R5):** done; shuffle did not collapse → "
            "neural typed aggregator **not** promoted.",
            "- **Age-primary seed-mask screen:** fold-0 panel audit **green** "
            "(`ok_for_seed_mask_gpu`); CUDA blocked only on 16-ep unlock — "
            "[`milestone-7g-prime-age-seed-mask.md`]"
            "(../../../docs/plans/milestone-7g-prime-age-seed-mask.md).",
            "- **Atlas association catalog:** done (SQL 013); non-blocking.",
            "- **Stage B CpG-panel GPU:** blocked until seed-mask screen + typed-RBS "
            "diagnostics.",
            "",
            "## Next",
            "",
            "Ordered ops:",
            "",
            "1. **Finish remaining 16-ep cascade queue** (scalar max→mean folds "
            "1–2 → vector mean→max ×3); do not kill GPU 0 jobs.",
            "2. **Let post-hoc enet finish** then re-sync `per_arm/` + "
            "`write_7g_gene_only_probe_report.py` / `apply_7g_16ep_decision.py`.",
            "3. **Age-primary seed-mask GPU** (`scripts/run_7g_prime_seed_mask.py "
            "--device cuda --reuse-panels`) — only after `promotion_decision.json` "
            "unlocks; fold 0, two seeds, K=256; do **not** launch Stage B.",
            "4. **Stage B** fold-selected CpG panel only after seed-mask screen.",
            "5. **Milestone 7** 5×6 OOF after Stage B + `direct_cpg.zarr`.",
            "",
        ]
    )
    (report_dir / "analysis.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--report-dir",
        type=Path,
        default=DEFAULT_REPORT,
        help="Inspection report directory",
    )
    args = parser.parse_args()
    paths = DataPaths.from_environment()
    report_dir = args.report_dir
    if not report_dir.is_absolute():
        report_dir = paths.project_root / report_dir
    report_dir.mkdir(parents=True, exist_ok=True)
    arms = load_per_arm(report_dir)
    cascade_rows: list[dict[str, Any]] = []
    cascade_folds_by_arm: dict[str, list[dict[str, Any]]] = {}
    for arm_id in CASCADE_ARMS:
        payload = arms.get(arm_id)
        if not payload:
            continue
        folds = payload.get("folds") or []
        cascade_folds_by_arm[arm_id] = folds
        e2e, _ = _cascade_tissue_f1(folds, "mbs_e2e")
        cascade_rows.append({"arm_id": arm_id, "mbs_e2e_f1": e2e})
    classical_rows: list[dict[str, Any]] = []
    cp = arms.get("C-mvalue-classical-G")
    if cp:
        for arm_name in CLASSICAL_SUFFIXES:
            f1, _ = _classical_tissue_f1(cp, arm_name)
            classical_rows.append({"arm_id": arm_name, "tissue_f1": f1})
    lock = build_lock_recommendation(
        cascade_rows,
        classical_rows,
        cascade_folds_by_arm=cascade_folds_by_arm,
        classical_payload=cp,
    )
    write_analysis(report_dir, lock=lock, paths=paths)
    print(f"wrote {report_dir / 'analysis.md'}", flush=True)


if __name__ == "__main__":
    main()
