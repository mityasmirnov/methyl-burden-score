#!/usr/bin/env python3
"""Milestone 12 OOF completeness gate (addresses a concrete review finding,
2026-09-08): the 5x6 runners (`run_12_nlight_oof.py`, `run_12_cascade_oof.py`)
treat any existing ``metrics.json`` as a completed run, and the nested-enet
post-hoc subprocess uses ``check=False`` -- a failed/NaN nested-enet fit can
silently coexist with a "completed" neural run. Neither runner's own summary
validates that all N=folds*restarts jobs have *valid* (present AND finite)
primary metrics before calling the campaign done.

This script is the missing gate: run it against a run-prefix before treating
a 5x6 as finished. Does not retrain anything -- read-only over existing
``metrics.json`` files.

Also flags (does not fix) a real mislabeling: ``primary_evaluation:
mbs_enet_nested`` in both OOF configs describes the *reported* readout, not
checkpoint *selection* -- neural checkpoint selection still uses
``validation_tissue_macro_f1_then_age_mae``. Nested enet is fit post-hoc on
whichever checkpoint that criterion picked; it cannot influence which epoch's
weights were kept.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mbs.paths import DataPaths  # noqa: E402

N_RESTARTS = 6


def _is_finite_num(x: Any) -> bool:
    return isinstance(x, int | float) and math.isfinite(x)


def _extract_metric_triplet(block: dict[str, Any] | None) -> dict[str, Any] | None:
    if not isinstance(block, dict):
        return None
    m = block.get("metrics")
    if not isinstance(m, dict):
        return None
    tissue_f1 = (m.get("tissue") or {}).get("macro_f1")
    age_mae = (m.get("age") or {}).get("mae")
    sex_auroc = (m.get("sex") or {}).get("auroc")
    return {"tissue_f1": tissue_f1, "age_mae": age_mae, "sex_auroc": sex_auroc}


def check_run(metrics_path: Path) -> dict[str, Any]:
    if not metrics_path.is_file():
        return {"status": "missing"}
    try:
        metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError) as exc:
        return {"status": "unreadable", "error": str(exc)}
    evaluations = metrics.get("evaluations") or {}
    e2e = _extract_metric_triplet(evaluations.get("mbs_e2e"))
    nested = _extract_metric_triplet(evaluations.get("mbs_enet_nested"))
    problems = []
    if e2e is None:
        problems.append("missing mbs_e2e")
    elif not all(_is_finite_num(v) for v in e2e.values()):
        problems.append(f"non-finite mbs_e2e: {e2e}")
    if nested is None:
        problems.append("missing mbs_enet_nested")
    elif not all(_is_finite_num(v) for v in nested.values()):
        problems.append(f"non-finite mbs_enet_nested: {nested}")
    if problems:
        return {"status": "invalid", "problems": problems, "e2e": e2e, "nested": nested}
    return {"status": "ok", "e2e": e2e, "nested": nested}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-prefix", required=True)
    parser.add_argument("--n-folds", type=int, required=True)
    parser.add_argument(
        "--fold-subdir",
        default=None,
        help="Cascade runs nest metrics under runs/<run_id>/fold_<i>/metrics.json; "
        "N-light writes runs/<run_id>/metrics.json directly. Pass 'fold' for the "
        "cascade layout, leave unset for the flat/N-light layout.",
    )
    args = parser.parse_args()

    paths = DataPaths.from_environment()
    n_jobs = args.n_folds * N_RESTARTS
    results: dict[str, dict[str, Any]] = {}
    for fold_i in range(args.n_folds):
        for restart in range(N_RESTARTS):
            run_id = f"{args.run_prefix}-f{fold_i}-r{restart}"
            run_root = paths.artifact_root / "runs" / run_id
            metrics_path = (
                run_root / f"fold_{fold_i}" / "metrics.json"
                if args.fold_subdir == "fold"
                else run_root / "metrics.json"
            )
            results[run_id] = check_run(metrics_path)

    n_ok = sum(1 for r in results.values() if r["status"] == "ok")
    print(f"[oof-gate] {args.run_prefix}: {n_ok}/{n_jobs} valid (fold x restart) results")
    for run_id, r in results.items():
        if r["status"] != "ok":
            print(f"[oof-gate]   {run_id}: {r['status']} {r.get('problems', r.get('error', ''))}")

    if n_ok < n_jobs:
        print(
            f"[oof-gate] INCOMPLETE: {n_jobs - n_ok} of {n_jobs} jobs are missing/invalid "
            "-- do not report this as a finished 5x6.",
            flush=True,
        )
        raise SystemExit(1)
    print("[oof-gate] COMPLETE: all jobs present with finite mbs_e2e and mbs_enet_nested.")


if __name__ == "__main__":
    main()
