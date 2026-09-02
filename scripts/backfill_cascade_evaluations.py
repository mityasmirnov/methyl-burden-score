#!/usr/bin/env python3
"""Backfill evaluations.mbs_e2e + fusion_full on pre-7G′ cascade runs (no retrain)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from mbs.paths import DataPaths
from mbs.training.cascade_loop import _evaluations_incomplete, train_cascade
from mbs.training.loop import load_experiment_config

ROOT = Path(__file__).resolve().parents[1]

BACKFILL_RUNS: tuple[tuple[str, str], ...] = (
    ("stage0-7g-cascade-v1", "configs/experiment/stage0_7g_methylation_eval.yaml"),
    ("stage0-7g-tissue-probe-P2", "configs/experiment/stage0_7g_cascade_tissue_probe_p2.yaml"),
    ("stage0-7g-tissue-probe-P4", "configs/experiment/stage0_7g_cascade_tissue_probe_p4.yaml"),
)


def _needs_backfill(paths: DataPaths, run_id: str) -> bool:
    run_root = paths.artifact_root / "runs" / run_id
    for fold_i in range(3):
        metrics_path = run_root / f"fold_{fold_i}" / "metrics.json"
        if not metrics_path.is_file():
            return True
        import json

        blob = json.loads(metrics_path.read_text(encoding="utf-8"))
        if _evaluations_incomplete(blob):
            return True
    return False


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--run-id", action="append", default=[])
    args = parser.parse_args()

    paths = DataPaths.from_environment()
    targets = BACKFILL_RUNS
    if args.run_id:
        targets = tuple((rid, "configs/experiment/stage0_7g_methylation_eval.yaml") for rid in args.run_id)

    for run_id, cfg_rel in targets:
        if not _needs_backfill(paths, run_id):
            print(f"[backfill] skip complete {run_id}", flush=True)
            continue
        cfg_path = paths.project_root / cfg_rel
        cfg = load_experiment_config(cfg_path)
        print(f"[backfill] eval_only {run_id} config={cfg_rel}", flush=True)
        train_cascade(
            project_root=paths.project_root,
            data_root=paths.data_root,
            artifact_root=paths.artifact_root,
            config=cfg,
            config_path=cfg_path,
            run_id=run_id,
            device_str=args.device,
            skip_if_done=True,
            eval_only=True,
        )
    print("[backfill] done", flush=True)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(130)
