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
from mbs.paths import DataPaths

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REPORT = ROOT / "reports/inspection/stage0_7g_gene_only_probe"

CASCADE_ARMS = ("P2-G", "P4-G", "P5-G-max", "P5-G-mean")
CLASSICAL_SUFFIXES = ("C-mvalue-ridge-G", "C-mvalue-enet-G", "C-mvalue-hgb-G", "C-mvalue-sva-G")
ARM_POOLING = {
    "P2-G": ("max", "max", 15),
    "P4-G": ("mean", "mean", 15),
    "P5-G-max": ("max", "max", 30),
    "P5-G-mean": ("mean", "mean", 30),
}
CLEAR_AHEAD_DELTA = 0.03


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


def _metric_from_fold(blob: dict[str, Any], metric_path: str) -> float | None:
    parts = metric_path.split(".")
    eval_keys = ("mbs_e2e", "mbs_linear_probe", "fusion_full", "fusion_mbs_direct")
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


def _cascade_tissue_f1(folds: list[dict[str, Any]], mode: str = "mbs_e2e") -> tuple[float | None, float | None]:
    path = f"{mode}.metrics.tissue.macro_f1"
    per = [_metric_from_fold(f, path) for f in folds]
    return _mean_std(per)


def _classical_tissue_f1(payload: dict[str, Any], arm_name: str) -> tuple[float | None, float | None]:
    per: list[float | None] = []
    for fold in payload.get("folds") or []:
        blob = (fold.get("arms") or {}).get(arm_name) or {}
        tissue = blob.get("tissue") or {}
        per.append(tissue.get("macro_f1"))
    return _mean_std(per)


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
) -> dict[str, Any]:
    best_cascade = max(cascade_rows, key=lambda r: r.get("mbs_e2e_f1") or -1.0, default=None)
    best_classical = max(classical_rows, key=lambda r: r.get("tissue_f1") or -1.0, default=None)
    rec: dict[str, Any] = {
        "primary_metric": "mbs_e2e.metrics.tissue.macro_f1",
        "locked_cascade_arm": None,
        "pooling_cpg": None,
        "pooling_region": None,
        "max_epochs": None,
        "cascade_clearly_ahead": None,
        "recommend_encoder_parity": False,
        "best_classical_arm": None,
    }
    if best_cascade and best_cascade.get("mbs_e2e_f1") is not None:
        arm_id = str(best_cascade["arm_id"])
        rec["locked_cascade_arm"] = arm_id
        pool = ARM_POOLING.get(arm_id)
        if pool:
            rec["pooling_cpg"], rec["pooling_region"], rec["max_epochs"] = pool
    if best_classical and best_classical.get("tissue_f1") is not None:
        rec["best_classical_arm"] = best_classical["arm_id"]
    c_f1 = (best_cascade or {}).get("mbs_e2e_f1")
    cl_f1 = (best_classical or {}).get("tissue_f1")
    if c_f1 is not None and cl_f1 is not None:
        rec["cascade_clearly_ahead"] = c_f1 >= cl_f1 + CLEAR_AHEAD_DELTA
        rec["recommend_encoder_parity"] = not rec["cascade_clearly_ahead"]
    return rec


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
        f"| Mode | Mean tissue macro-F1 |",
        f"|------|---------------------:|",
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


def write_analysis(report_dir: Path, *, lock: dict[str, Any]) -> None:
    arms = load_per_arm(report_dir)
    cascade_rows: list[dict[str, Any]] = []
    for arm_id in CASCADE_ARMS:
        payload = arms.get(arm_id)
        if payload is None:
            continue
        folds = payload.get("folds") or []
        e2e, e2e_std = _cascade_tissue_f1(folds, "mbs_e2e")
        probe, _ = _cascade_tissue_f1(folds, "mbs_linear_probe")
        cascade_rows.append(
            {
                "arm_id": arm_id,
                "mbs_e2e_f1": e2e,
                "mbs_e2e_std": e2e_std,
                "mbs_linear_probe_f1": probe,
                "n_folds": len(folds),
            }
        )

    classical_payload = arms.get("C-mvalue-classical-G")
    classical_rows: list[dict[str, Any]] = []
    if classical_payload:
        for arm_name in CLASSICAL_SUFFIXES:
            f1, f1_std = _classical_tissue_f1(classical_payload, arm_name)
            classical_rows.append(
                {"arm_id": arm_name, "tissue_f1": f1, "tissue_f1_std": f1_std}
            )

    if not lock:
        lock = build_lock_recommendation(cascade_rows, classical_rows)
    write_json(report_dir / "lock_recommendation.json", lock)

    lines = [
        "# 7G′ Stage A — gene-only MBS architecture selection",
        "",
        "Primary metric: **`mbs_e2e`** tissue macro-F1 (end-to-end MBS heads; not late fusion).",
        "Classical comparator: **`C-mvalue-*-G`** on identical `gene_cols`.",
        "",
        "## Cascade arms (gene-linked CpGs only)",
        "",
        "| Arm | mbs_e2e F1 | mbs_linear_probe F1 | folds |",
        "|-----|-----------:|--------------------:|------:|",
    ]
    for row in sorted(cascade_rows, key=lambda r: r.get("mbs_e2e_f1") or -1.0, reverse=True):
        lines.append(
            f"| {row['arm_id']} | {_fmt(row.get('mbs_e2e_f1'))} "
            f"(±{_fmt(row.get('mbs_e2e_std'))}) | "
            f"{_fmt(row.get('mbs_linear_probe_f1'))} | {row.get('n_folds', 0)} |"
        )
    lines.extend(
        [
            "",
            "## Classical arms (-G panel)",
            "",
            "| Arm | tissue macro-F1 |",
            "|-----|----------------:|",
        ]
    )
    for row in sorted(classical_rows, key=lambda r: r.get("tissue_f1") or -1.0, reverse=True):
        lines.append(f"| {row['arm_id']} | {_fmt(row.get('tissue_f1'))} (±{_fmt(row.get('tissue_f1_std'))}) |")

    lines.extend(
        [
            "",
            "## Locked architecture (Stage B input)",
            "",
            f"- **Cascade arm:** `{lock.get('locked_cascade_arm')}`",
            f"- **Pooling (CpG / region):** `{lock.get('pooling_cpg')}` / `{lock.get('pooling_region')}`",
            f"- **Epoch ceiling:** {lock.get('max_epochs')}",
            f"- **Best classical:** `{lock.get('best_classical_arm')}`",
            f"- **Cascade clearly ahead (≥{CLEAR_AHEAD_DELTA} F1):** {lock.get('cascade_clearly_ahead')}",
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

    lines.append(orphan_ablation_section(arms.get("P2-orphan-ablation")))
    lines.extend(
        [
            "## Next",
            "",
            "- Stage B: fold-safe `C-mvalue-enetS`, `N-cascade-S`, `N-light-type` (FlatDeepSetRegion), "
            "`N-full` / `N-mbs-direct-only`, plus `direct_cpg.zarr`.",
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
    lock = build_lock_recommendation([], [])
    arms = load_per_arm(report_dir)
    cascade_rows: list[dict[str, Any]] = []
    for arm_id in CASCADE_ARMS:
        payload = arms.get(arm_id)
        if not payload:
            continue
        e2e, _ = _cascade_tissue_f1(payload.get("folds") or [], "mbs_e2e")
        cascade_rows.append({"arm_id": arm_id, "mbs_e2e_f1": e2e})
    classical_rows: list[dict[str, Any]] = []
    cp = arms.get("C-mvalue-classical-G")
    if cp:
        for arm_name in CLASSICAL_SUFFIXES:
            f1, _ = _classical_tissue_f1(cp, arm_name)
            classical_rows.append({"arm_id": arm_name, "tissue_f1": f1})
    lock = build_lock_recommendation(cascade_rows, classical_rows)
    write_analysis(report_dir, lock=lock)
    print(f"wrote {report_dir / 'analysis.md'}", flush=True)


if __name__ == "__main__":
    main()
