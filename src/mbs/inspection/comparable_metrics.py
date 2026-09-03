"""Panel × eval-mode tissue F1 rows for fair 7G / 7G′ comparisons."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

# Each spec: arm label, run_id under $MBS_ARTIFACT_ROOT/runs/, panel, eval_mode, dotted metric path.
COMPARISON_SPECS: tuple[dict[str, str], ...] = (
    {
        "arm": "P0-baseline",
        "run_id": "stage0-7g-cascade-v1",
        "panel": "65k prefix",
        "eval_mode": "fusion_full",
        "metric_path": "metrics.tissue.macro_f1",
    },
    {
        "arm": "P0-baseline",
        "run_id": "stage0-7g-cascade-v1",
        "panel": "65k prefix",
        "eval_mode": "mbs_e2e",
        "metric_path": "evaluations.mbs_e2e.metrics.tissue.macro_f1",
    },
    {
        "arm": "P2-end2end-tissue-weight",
        "run_id": "stage0-7g-tissue-probe-P2",
        "panel": "65k prefix",
        "eval_mode": "fusion_full",
        "metric_path": "metrics.tissue.macro_f1",
    },
    {
        "arm": "P2-end2end-tissue-weight",
        "run_id": "stage0-7g-tissue-probe-P2",
        "panel": "65k prefix",
        "eval_mode": "mbs_e2e",
        "metric_path": "evaluations.mbs_e2e.metrics.tissue.macro_f1",
    },
    {
        "arm": "P4-pooling-mean",
        "run_id": "stage0-7g-tissue-probe-P4",
        "panel": "65k prefix",
        "eval_mode": "fusion_full",
        "metric_path": "metrics.tissue.macro_f1",
    },
    {
        "arm": "P4-pooling-mean",
        "run_id": "stage0-7g-tissue-probe-P4",
        "panel": "65k prefix",
        "eval_mode": "mbs_e2e",
        "metric_path": "evaluations.mbs_e2e.metrics.tissue.macro_f1",
    },
    {
        "arm": "P5-epochs-30",
        "run_id": "stage0-7g-tissue-probe-P5",
        "panel": "65k prefix",
        "eval_mode": "fusion_full",
        "metric_path": "metrics.tissue.macro_f1",
    },
    {
        "arm": "P5-epochs-30",
        "run_id": "stage0-7g-tissue-probe-P5",
        "panel": "65k prefix",
        "eval_mode": "mbs_e2e",
        "metric_path": "evaluations.mbs_e2e.metrics.tissue.macro_f1",
    },
    {
        "arm": "P2-G",
        "run_id": "stage0-7g-gene-probe-P2-G-explicit",
        "panel": "gene-linked",
        "eval_mode": "fusion_full",
        "metric_path": "evaluations.fusion_full.metrics.tissue.macro_f1",
    },
    {
        "arm": "P2-G",
        "run_id": "stage0-7g-gene-probe-P2-G-explicit",
        "panel": "gene-linked",
        "eval_mode": "mbs_e2e",
        "metric_path": "evaluations.mbs_e2e.metrics.tissue.macro_f1",
    },
    {
        "arm": "P2-G",
        "run_id": "stage0-7g-gene-probe-P2-G-explicit",
        "panel": "gene-linked",
        "eval_mode": "mbs_enet",
        "metric_path": "evaluations.mbs_enet.metrics.tissue.macro_f1",
    },
    {
        "arm": "P4-G",
        "run_id": "stage0-7g-gene-probe-P4-G-explicit",
        "panel": "gene-linked",
        "eval_mode": "fusion_full",
        "metric_path": "evaluations.fusion_full.metrics.tissue.macro_f1",
    },
    {
        "arm": "P4-G",
        "run_id": "stage0-7g-gene-probe-P4-G-explicit",
        "panel": "gene-linked",
        "eval_mode": "mbs_e2e",
        "metric_path": "evaluations.mbs_e2e.metrics.tissue.macro_f1",
    },
    {
        "arm": "P4-G",
        "run_id": "stage0-7g-gene-probe-P4-G-explicit",
        "panel": "gene-linked",
        "eval_mode": "mbs_enet",
        "metric_path": "evaluations.mbs_enet.metrics.tissue.macro_f1",
    },
    {
        "arm": "P5-G-max",
        "run_id": "stage0-7g-gene-probe-P5-G-max-explicit",
        "panel": "gene-linked",
        "eval_mode": "mbs_e2e",
        "metric_path": "evaluations.mbs_e2e.metrics.tissue.macro_f1",
    },
    {
        "arm": "C-mvalue-enet",
        "run_id": "",
        "panel": "65k prefix",
        "eval_mode": "classical",
        "metric_path": "classical:C-mvalue-enet",
    },
    {
        "arm": "C-mvalue-enet-G",
        "run_id": "",
        "panel": "gene-linked",
        "eval_mode": "classical",
        "metric_path": "classical:C-mvalue-enet-G",
    },
)


def _metric_from_blob(blob: dict[str, Any], metric_path: str) -> float | None:
    parts = metric_path.split(".")
    cur: Any = blob
    for key in parts:
        if not isinstance(cur, dict):
            return None
        cur = cur.get(key)
    if cur is None or (isinstance(cur, float) and math.isnan(cur)):
        return None
    return float(cur)


def _mean_std(vals: list[float | None]) -> tuple[float | None, float | None]:
    nums = [float(v) for v in vals if v is not None and not math.isnan(v)]
    if not nums:
        return None, None
    arr = np.asarray(nums, dtype=np.float64)
    return float(arr.mean()), float(arr.std(ddof=1)) if arr.size > 1 else 0.0


def _classical_f1_per_fold(
    classical_path: Path,
    arm_name: str,
) -> list[float | None]:
    if not classical_path.is_file():
        return []
    payload = json.loads(classical_path.read_text(encoding="utf-8"))
    out: list[float | None] = []
    for fold in payload.get("folds") or []:
        blob = (fold.get("arms") or {}).get(arm_name) or {}
        tissue = blob.get("tissue") if isinstance(blob.get("tissue"), dict) else {}
        f1 = tissue.get("macro_f1")
        out.append(float(f1) if f1 is not None else None)
    return out


def _mbs_e2e_fold_valid(blob: dict[str, Any]) -> bool:
    evaluations = blob.get("evaluations") or {}
    e2e = evaluations.get("mbs_e2e")
    if not isinstance(e2e, dict):
        return False
    return e2e.get("eval_split") == "test"


def load_comparable_rows(
    artifact_root: Path,
    *,
    classical_baselines_path: Path | None = None,
    gene_classical_path: Path | None = None,
) -> list[dict[str, Any]]:
    """Aggregate tissue macro-F1 per arm × panel × eval_mode."""
    rows: list[dict[str, Any]] = []
    for spec in COMPARISON_SPECS:
        per_fold: list[float | None] = []
        if spec["metric_path"].startswith("classical:"):
            arm_name = spec["metric_path"].split(":", 1)[1]
            classical_path = (
                gene_classical_path
                if arm_name.endswith("-G") and gene_classical_path is not None
                else classical_baselines_path
            )
            if classical_path is not None:
                per_fold = _classical_f1_per_fold(classical_path, arm_name)
        elif spec["run_id"]:
            run_root = artifact_root / "runs" / spec["run_id"]
            for fold_i in range(3):
                metrics_path = run_root / f"fold_{fold_i}" / "metrics.json"
                if not metrics_path.is_file():
                    per_fold.append(None)
                    continue
                blob = json.loads(metrics_path.read_text(encoding="utf-8"))
                if spec["eval_mode"] == "mbs_e2e" and not _mbs_e2e_fold_valid(blob):
                    per_fold.append(None)
                    continue
                per_fold.append(_metric_from_blob(blob, spec["metric_path"]))
        mean_f1, std_f1 = _mean_std(per_fold)
        row: dict[str, Any] = {
            "arm": spec["arm"],
            "run_id": spec["run_id"],
            "panel": spec["panel"],
            "eval_mode": spec["eval_mode"],
            "tissue_macro_f1": mean_f1,
            "tissue_macro_f1_std": std_f1,
            "n_folds": sum(1 for v in per_fold if v is not None),
            "per_fold_f1": per_fold,
        }
        if spec["eval_mode"] == "mbs_e2e" and row["n_folds"] == 0 and spec["run_id"]:
            row["invalid_reason"] = "mbs_e2e missing eval_split=test (train+test leak)"
        rows.append(row)
    return rows


def render_comparable_ranking_section(rows: list[dict[str, Any]]) -> list[str]:
    """Markdown: fair comparison table (only compare rows with same panel + eval_mode)."""
    lines = [
        "## Comparable ranking (panel × eval mode)",
        "",
        "Compare **only within the same row group** (same panel and eval mode). "
        "Stage A primary metric is **`mbs_e2e`** on the **gene-linked** panel (test split only). "
        "7G tissue probe P0–P5 used **late fusion (`fusion_full`)** on the **65k prefix**. "
        "Rows marked *invalid* used pre-fix `mbs_e2e` that scored train+validation+test together.",
        "",
        "| Arm | Panel | Eval mode | Tissue macro-F1 | folds | Notes |",
        "|-----|-------|-----------|----------------:|------:|-------|",
    ]
    for row in rows:
        f1 = row.get("tissue_macro_f1")
        std = row.get("tissue_macro_f1_std")
        if f1 is None:
            disp = "—"
        elif std is not None and row.get("n_folds", 0) > 1:
            disp = f"{f1:.3f} (±{std:.3f})"
        else:
            disp = f"{f1:.3f}"
        note = row.get("invalid_reason") or ""
        if note:
            disp = f"*{disp}*" if disp != "—" else "—"
        lines.append(
            f"| `{row['arm']}` | {row['panel']} | `{row['eval_mode']}` | {disp} | "
            f"{row.get('n_folds', 0)} | {note} |"
        )
    lines.extend(
        [
            "",
            "**Fair pairs (examples):**",
            "",
            "- **Late fusion, 65k:** `P2-end2end` vs `P4-pooling-mean` vs `C-mvalue-enet`.",
            "- **Late fusion, gene-linked:** `P2-G` / `P4-G` `fusion_full` on `explicit_only`.",
            "- **MBS e2e, gene-linked:** `P2-G` vs `P4-G` vs `P5-G-max` vs `C-mvalue-enet-G` (Stage A lock).",
            "- **MBS enet readout:** same frozen MBS as e2e, elastic-net heads (not a lock substitute).",
            "",
        ]
    )
    return lines
