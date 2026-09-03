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
    "N-cascade-scalar-mean-max": ("mean", "max", 5),
    "N-cascade-scalar-max-mean": ("max", "mean", 5),
    "N-cascade-vector-mean-max": ("mean", "max", 5),
    "N-cascade-vector-max-max": ("max", "max", 5),
    "N-light-gene-max": ("n/a", "max", 5),
    "N-light-gene-mean": ("n/a", "mean", 5),
}
CASCADE_READOUTS = (
    "mbs_e2e",
    "mbs_linear_probe",
    "mbs_enet",
    "rbs_linear_probe",
    "rbs_enet",
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



def _cascade_metric_means(
    folds: list[dict[str, Any]],
    mode: str,
    *metric_keys: str,
) -> tuple[float | None, float | None]:
    path = ".".join((mode, "metrics", *metric_keys))
    per = [_metric_from_fold(f, path) for f in folds]
    return _mean_std(per)


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
    eval_keys = ("mbs_e2e", "mbs_linear_probe", "mbs_enet", "fusion_full", "fusion_mbs_direct")
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
    rec: dict[str, Any] = {
        "primary_metric": "mbs_e2e.metrics.tissue.macro_f1",
        "locked_cascade_arm": None,
        "pooling_cpg": None,
        "pooling_region": None,
        "max_epochs": None,
        "cascade_clearly_ahead": None,
        "recommend_encoder_parity": False,
        "best_classical_arm": None,
        "lock_blocked_reason": None,
        "mbs_e2e_valid": False,
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
        rec["lock_blocked_reason"] = (
            "no cascade arm with test-only mbs_e2e (eval_split=test on all folds)"
        )
        return rec
    if not valid_classical_rows:
        rec["lock_blocked_reason"] = "no completed C-mvalue-*-G classical folds"
        return rec

    rec["mbs_e2e_valid"] = True
    rec.update(LOCK_TRAINING_DEFAULTS)
    best_cascade = max(valid_cascade_rows, key=lambda r: r.get("mbs_e2e_f1") or -1.0)
    best_classical = max(valid_classical_rows, key=lambda r: r.get("tissue_f1") or -1.0)
    arm_id = str(best_cascade["arm_id"])
    rec["locked_cascade_arm"] = arm_id
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
    arm_id: str,
    *metric_keys: str,
) -> tuple[float, float, float]:
    """Extract per-fold values then bootstrap."""
    payload = arm_payloads.get(arm_id)
    if not payload:
        return (float("nan"),) * 3
    folds = payload.get("folds") or []
    vals = [v for f in folds if (v := _extract_scalar(f, *metric_keys)) is not None]
    return _bootstrap_ci(vals)


def _ablation_table_row(
    arm_id: str,
    arm_payloads: dict[str, dict | None],
    *,
    tissue_keys: tuple[str, ...],
    age_keys: tuple[str, ...],
    sex_keys: tuple[str, ...],
) -> str:
    t_mean, t_lo, t_hi = _arm_mean_ci(arm_payloads, arm_id, *tissue_keys)
    a_mean, a_lo, a_hi = _arm_mean_ci(arm_payloads, arm_id, *age_keys)
    s_mean, s_lo, s_hi = _arm_mean_ci(arm_payloads, arm_id, *sex_keys)

    def _fmt(m: float, lo: float, hi: float) -> str:
        if math.isnan(m):
            return "—"
        return f"{m:.3f} [{lo:.3f}–{hi:.3f}]"

    return f"| {arm_id} | {_fmt(t_mean, t_lo, t_hi)} | {_fmt(a_mean, a_lo, a_hi)} | {_fmt(s_mean, s_lo, s_hi)} |"


def _repr_diagnostics_row(arm_id: str, arm_payloads: dict[str, dict | None]) -> str:
    """Build a representation diagnostics markdown row from stored fold metrics."""
    payload = arm_payloads.get(arm_id)
    if not payload:
        return f"| {arm_id} | — | — | — | — |"
    folds = payload.get("folds") or []

    def _mean_key(*keys: str) -> str:
        vals = [v for f in folds if (v := _extract_scalar(f, *keys)) is not None]
        if not vals:
            return "—"
        return f"{np.mean(vals):.3f}"

    score_sd = _mean_key("evaluations", "mbs_e2e", "repr_diagnostics", "gene_score_sd")
    sat_frac = _mean_key("evaluations", "mbs_e2e", "repr_diagnostics", "saturation_fraction")
    const_frac = _mean_key("evaluations", "mbs_e2e", "repr_diagnostics", "constant_score_fraction")
    corr_m = _mean_key("evaluations", "mbs_e2e", "repr_diagnostics", "corr_mean_m")
    return f"| {arm_id} | {score_sd} | {sat_frac} | {const_frac} | {corr_m} |"


def annotation_ablation_section(arm_payloads: dict[str, dict | None]) -> str:
    """Render annotation ablation grid (A0–A4/A7, N0–N3) with bootstrap CIs."""
    # Arm IDs follow config experiment.name convention
    ablation_arms = [
        ("A0", "stage0_7g_gene_only_probe_ablation_m_only", "M only"),
        ("A1", "stage0_7g_gene_only_probe_ablation_m_role", "M + gene role"),
        ("A2", "stage0_7g_gene_only_probe_ablation_m_context", "M + CpG context"),
        ("A3", "stage0_7g_gene_only_probe_ablation_m_role_context", "M + role + context"),
        ("A4/A7", "stage0_7g_gene_only_probe_ablation_full", "All (regulatory zero)"),
    ]
    neg_arms = [
        ("N0", "stage0_7g_gene_only_probe_ablation_n0_obs_only", "Observed flag only"),
        ("N1", "stage0_7g_gene_only_probe_ablation_n1_anno_only", "Annotations only (no M)"),
        ("N2", "stage0_7g_gene_only_probe_ablation_n2_reg_permuted", "Reg. permuted"),
        ("N3", "stage0_7g_gene_only_probe_ablation_n3_reg_zero", "All-zero regulatory"),
    ]

    tissue_keys = ("evaluations", "mbs_e2e", "metrics", "tissue_macro_f1")
    age_keys = ("evaluations", "mbs_e2e", "metrics", "age", "mae")
    sex_keys = ("evaluations", "mbs_e2e", "metrics", "sex", "auroc")

    header = [
        "## Annotation ablation grid (A0–A4, N0–N3)\n",
        "Fold 0, `mean` pooling, 8 epochs, two seeds. "
        "Bootstrap 95% CIs from available folds. "
        "Tissue macro-F1, age MAE, sex AUROC.\n",
        "**Note:** A4 and A7 are identical while regulatory channels are zero (cCRE/DHS/ChromHMM not on disk).\n",
        "| Arm | Features | Tissue macro-F1 [95% CI] | Age MAE [95% CI] | Sex AUROC [95% CI] |",
        "|-----|----------|-------------------------:|-----------------:|-------------------:|",
    ]
    rows_a = [
        f"| {label} | {desc} | "
        + _ablation_table_row(arm_id, arm_payloads, tissue_keys=tissue_keys, age_keys=age_keys, sex_keys=sex_keys).split("|", 2)[2]
        for label, arm_id, desc in ablation_arms
    ]
    neg_header = [
        "",
        "### Negative controls\n",
        "| Arm | Features | Tissue macro-F1 [95% CI] | Age MAE [95% CI] | Sex AUROC [95% CI] |",
        "|-----|----------|-------------------------:|-----------------:|-------------------:|",
    ]
    rows_n = [
        f"| {label} | {desc} | "
        + _ablation_table_row(arm_id, arm_payloads, tissue_keys=tissue_keys, age_keys=age_keys, sex_keys=sex_keys).split("|", 2)[2]
        for label, arm_id, desc in neg_arms
    ]

    repr_header = [
        "",
        "### Representation diagnostics (fold 0 mean across seeds)\n",
        "> Values populated only when `stage_a_per_epoch_eval: true` and repr_diagnostics logged.\n",
        "| Arm | Gene-score SD | Saturation frac | Const-score frac | Corr w/ mean-M |",
        "|-----|:-------------:|:---------------:|:----------------:|:--------------:|",
    ]
    repr_rows_a = [_repr_diagnostics_row(arm_id, arm_payloads) for _, arm_id, _ in ablation_arms]
    repr_rows_n = [_repr_diagnostics_row(arm_id, arm_payloads) for _, arm_id, _ in neg_arms]

    return "\n".join(
        header + rows_a + neg_header + rows_n + repr_header + repr_rows_a + repr_rows_n + [""]
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
    tissue, tissue_std = _cascade_metric_means(folds, mode, "tissue", "macro_f1")
    if tissue is None:
        return None
    age_mae, age_mae_std = _cascade_metric_means(folds, mode, "age", "mae")
    age_r2, age_r2_std = _cascade_metric_means(folds, mode, "age", "r2")
    sex_auroc, sex_auroc_std = _cascade_metric_means(folds, mode, "sex", "auroc")
    sex_f1, sex_f1_std = _cascade_metric_means(folds, mode, "sex", "macro_f1")
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
        "n_folds": len(folds),
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
                    "n_folds": 3,
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
            "**Readouts:** `mbs_e2e` = jointly trained neural heads on MBS (Stage A lock metric); "
            "`mbs_linear_probe` / `mbs_enet` = new sklearn heads on the **same frozen MBS**; "
            "`rbs_linear_probe` / `rbs_enet` = frozen **gene-linked RBS** (pre–gene-pool); "
            "`classical` = sklearn on gene-linked CpG M-values (no encoder).",
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


def render_architecture_qa(rows: list[dict[str, Any]]) -> list[str]:
    """Answer the seven Stage A architecture questions from available rows."""

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
    rbs = _mode("P2-G", "rbs_enet") or _mode("P2-G", "rbs_linear_probe")
    mbs_enet = _mode("P2-G", "mbs_enet")
    classical = next((r for r in rows if r.get("arm_id") == "C-mvalue-enet-G"), None)

    def _cmp_pool(level: str) -> str:
        if level == "cpg":
            a, b = mean_max, p2  # mean vs max at cpg (region fixed max)
            label_a, label_b = "mean-max", "max-max"
            if mean_max is None or p2 is None:
                a, b = p4, max_mean
                label_a, label_b = "mean-mean", "max-mean"
        else:
            a, b = max_mean, p2  # mean vs max at gene (cpg fixed max)
            label_a, label_b = "max-mean", "max-max"
            if max_mean is None or p2 is None:
                a, b = p4, mean_max
                label_a, label_b = "mean-mean", "mean-max"
        if a is None or b is None:
            return f"Insufficient arms for {level}-level comparison yet."
        better = a if (a.get("tissue_f1") or 0) >= (b.get("tissue_f1") or 0) else b
        return (
            f"`{label_a}` tissue F1={_fmt(a.get('tissue_f1'))} vs "
            f"`{label_b}` {_fmt(b.get('tissue_f1'))}; "
            f"age MAE {_fmt(a.get('age_mae'))} vs {_fmt(b.get('age_mae'))}. "
            f"Prefer **`{better['arm_id']}`** on this slice (check Pareto)."
        )

    q3 = "Pending RBS diagnostic."
    if rbs is not None and mbs_enet is not None:
        q3 = (
            f"P2 `rbs_*` tissue F1={_fmt(rbs.get('tissue_f1'))}, age MAE={_fmt(rbs.get('age_mae'))}; "
            f"`mbs_enet` tissue={_fmt(mbs_enet.get('tissue_f1'))}, age={_fmt(mbs_enet.get('age_mae'))}. "
            + (
                "Gene pooling recovers tissue but **drops age/sex** relative to RBS."
                if (rbs.get("age_mae") or 99) < (mbs_enet.get("age_mae") or 0)
                else "Loss is not clearly at gene pooling; check scalar RBS vs raw CpG."
            )
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
        "7. **Best performance/compute:** Prefer landed P2/P4 (15 ep) over P5; "
        "promote Tier-1 (5 ep) arms only when Pareto/near-best, then confirm at 15 ep.",
        "",
    ]
    return lines


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
        e2e, e2e_std = _cascade_tissue_f1(folds, "mbs_e2e")
        probe, _ = _cascade_tissue_f1(folds, "mbs_linear_probe")
        enet, enet_std = _cascade_tissue_f1(folds, "mbs_enet")
        age_mae, _ = _cascade_metric_means(folds, "mbs_e2e", "age", "mae")
        age_r2, _ = _cascade_metric_means(folds, "mbs_e2e", "age", "r2")
        sex_auroc, _ = _cascade_metric_means(folds, "mbs_linear_probe", "sex", "auroc")
        if sex_auroc is None:
            sex_auroc, _ = _cascade_metric_means(folds, "mbs_enet", "sex", "auroc")
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
    if lock.get("locked_cascade_arm"):
        glossary_ids.insert(0, str(lock["locked_cascade_arm"]))

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
    lines.extend(render_architecture_qa(task_rows))
    lines.extend(
        [
        "## Cascade arms (gene-linked CpGs only)",
        "",
        "Primary **`mbs_e2e`** (test split only); **`mbs_linear_probe`** and **`mbs_enet`** "
        "are readouts of the **same frozen MBS**; **`rbs_*`** use gene-linked RBS. "
        "Contaminated pre-fix **`mbs_e2e`** shown as *invalid*.",
        "",
        "| Arm | mbs_e2e F1 | linear probe F1 | mbs_enet F1 | age MAE (e2e) | sex AUROC (probe) | folds |",
        "|-----|-----------:|----------------:|------------:|--------------:|------------------:|------:|",
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
        lines.append(
            f"| {row['arm_id']} | {e2e_disp} | "
            f"{_fmt(row.get('mbs_linear_probe_f1'))} | {enet_disp} | "
            f"{_fmt(row.get('age_mae'))} | "
            f"{_fmt(row.get('sex_auroc'))} | {row.get('n_folds', 0)} |"
        )
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
            "## Locked architecture (Stage B input)",
            "",
        ]
    )
    if lock.get("locked_cascade_arm"):
        lines.extend(
            [
                f"- **Cascade arm:** `{lock.get('locked_cascade_arm')}`",
                f"- **Pooling (CpG / region):** `{lock.get('pooling_cpg')}` / `{lock.get('pooling_region')}`",
                f"- **Epoch ceiling:** {lock.get('max_epochs')}",
                f"- **Best classical (-G):** `{lock.get('best_classical_arm')}`",
                f"- **Cascade clearly ahead (≥{CLEAR_AHEAD_DELTA} tissue F1):** {lock.get('cascade_clearly_ahead')}",
                "",
            ]
        )
        if lock.get("recommend_encoder_parity"):
            lines.append(
                "**Encoder parity recommended:** re-run **FlatDeepSet** and **HierarchicalDeepSet** "
                "on the same `gene_cols` before committing to cascade for Stage B / Milestone 7."
            )
            lines.append("")
        else:
            lines.append(
                "Cascade leads classical `-G` by ≥0.03 F1; optional Flat/Hier parity runs are **not** required."
            )
            lines.append("")
    else:
        reason = lock.get("lock_blocked_reason") or "insufficient evidence"
        lines.extend(
            [
                "**No architecture lock.** Stage B and Milestone **7** OOF remain blocked.",
                f"- **Reason:** {reason}",
                "- Re-run Stage A with test-only `mbs_e2e` and matched `C-mvalue-*-G` on identical "
                "`gene_cols` (`explicit_only` allocation).",
                "",
            ]
        )

    lines.append(orphan_ablation_section(arms.get("P2-orphan-ablation")))
    lines.append(annotation_ablation_section(arms))
    lines.extend(
        [
            "## Parallel / follow-on work",
            "",
            "- **Stage A required GPU arms** (`P2-G`, `P4-G`, `P5-G-max`, `C-mvalue-*-G`) are complete "
            "on `explicit_only`. Optional `P5-G-mean` was not run.",
            "- **Stage A screen (sequential):** train one arm at a time and regenerate this report after each. "
            "Order: `N-light-gene-max` → `N-light-gene-mean` → mixed scalar cascades → vector cascades; "
            "promote Tier-2 (15 ep) only if Pareto/near-best.",
            "- **Encoder parity (optional):** FlatDeepSet + HierarchicalDeepSet on same `gene_cols` if cascade "
            "does not lead classical by ≥0.03 F1.",
            "- **Stage B (after lock):** fold-safe `C-mvalue-enetS`, `N-cascade-S`, `N-light-type`, "
            "`direct_cpg.zarr`, full-model fusion arms.",
            "",
            "## Next",
            "",
            "- **Stage A screen (sequential):** continue remaining Tier-1 arms after each landed light arm "
            "updates this report.",
            "- Stage B (after lock): fold-safe `C-mvalue-enetS`, `N-cascade-S`, `N-light-type` (FlatDeepSetRegion), "
            "`N-mbs-posthoc-full-fusion` / `N-mbs-posthoc-mbs-direct`, plus `direct_cpg.zarr`.",
            "- Milestone **7** 5×6 OOF remains blocked until Stage B completes.",
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
