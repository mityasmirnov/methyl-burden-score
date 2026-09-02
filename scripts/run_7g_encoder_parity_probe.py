#!/usr/bin/env python3
"""Optional 7G′ encoder parity: FlatDeepSet + HierarchicalDeepSet on gene_cols."""

from __future__ import annotations

import argparse
import json
import sys
from copy import deepcopy
from pathlib import Path
from typing import Any

import numpy as np

from mbs.annotation.manifest import write_json
from mbs.paths import DataPaths
from mbs.training.dev_cv import inject_fold_into_config, load_frozen_folds
from mbs.training.hier_loop import train_hierarchical_baseline
from mbs.training.loop import load_experiment_config, train_flat_baseline

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "configs/experiment/stage0_7g_encoder_parity.yaml"


def _fold_run_id(prefix: str, fold_idx: int) -> str:
    return f"{prefix}-f{fold_idx}"


def _metrics_done(run_root: Path) -> bool:
    return (run_root / "metrics.json").is_file()


def _tissue_f1(metrics: dict[str, Any]) -> float | None:
    tissue = metrics.get("tissue")
    if isinstance(tissue, dict) and tissue.get("macro_f1") is not None:
        return float(tissue["macro_f1"])
    if metrics.get("macro_f1") is not None:
        return float(metrics["macro_f1"])
    return None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--force", action="store_true", help="Run even if lock says skip")
    parser.add_argument("--skip-if-done", action="store_true", default=True)
    parser.add_argument("--no-skip-if-done", dest="skip_if_done", action="store_false")
    parser.add_argument("--arm", action="append", default=[])
    args = parser.parse_args()

    paths = DataPaths.from_environment()
    config_path = args.config if args.config.is_absolute() else paths.project_root / args.config
    cfg = load_experiment_config(config_path)
    report_rel = Path(str(cfg.get("report_dir", "reports/inspection/stage0_7g_encoder_parity")))
    report_dir = report_rel if report_rel.is_absolute() else paths.project_root / report_rel
    report_dir.mkdir(parents=True, exist_ok=True)

    lock_path = paths.project_root / str(
        cfg.get(
            "lock_from_stage_a",
            "reports/inspection/stage0_7g_gene_only_probe/lock_recommendation.json",
        )
    )
    lock: dict[str, Any] = {}
    if lock_path.is_file():
        lock = json.loads(lock_path.read_text(encoding="utf-8"))
    if lock.get("recommend_encoder_parity") is False and not args.force:
        print(
            "[encoder-parity] skip: lock recommends against parity "
            f"(cascade leads classical ≥0.03 F1); use --force to override",
            flush=True,
        )
        return

    split_id = str(cfg.get("split_id", "hub-ats-7e-3fold-v1"))
    max_loci = int(cfg.get("cv_budget", {}).get("max_loci", 65536))
    max_epochs = int(cfg.get("cv_budget", {}).get("max_epochs", 15))
    folds_path = paths.artifact_root / "splits" / split_id / "folds.json"
    fold_pack = load_frozen_folds(folds_path)

    arms = cfg.get("arms") or []
    requested = set(args.arm)
    results: dict[str, Any] = {
        "milestone": "7G-encoder-parity",
        "lock_from_stage_a": lock,
        "arms": {},
    }

    for arm in arms:
        arm_id = str(arm["id"])
        if requested and arm_id not in requested:
            continue
        kind = str(arm["kind"])
        run_prefix = str(arm["run_prefix"])
        rel_cfg = Path(str(arm["config"]))
        arm_cfg_path = rel_cfg if rel_cfg.is_absolute() else paths.project_root / rel_cfg
        base_cfg = load_experiment_config(arm_cfg_path)
        fold_rows: list[dict[str, Any]] = []

        for fold_idx, fold in enumerate(fold_pack["folds"]):
            run_id = _fold_run_id(run_prefix, fold_idx)
            run_root = paths.artifact_root / "runs" / run_id
            if args.skip_if_done and _metrics_done(run_root):
                metrics = json.loads((run_root / "metrics.json").read_text(encoding="utf-8"))
                fold_rows.append({"fold": fold_idx, "run_id": run_id, "metrics": metrics, "resumed": True})
                print(f"[encoder-parity] skip done {arm_id} fold={fold_idx}", flush=True)
                continue

            fold_cfg = inject_fold_into_config(deepcopy(base_cfg), fold, seed=42 + fold_idx)
            fold_cfg.setdefault("cv_budget", {})["max_loci"] = max_loci
            fold_cfg.setdefault("training", {})["max_epochs"] = max_epochs
            print(f"[encoder-parity] train {arm_id} fold={fold_idx} run_id={run_id}", flush=True)
            if kind == "flat_train":
                result = train_flat_baseline(
                    project_root=paths.project_root,
                    data_root=paths.data_root,
                    artifact_root=paths.artifact_root,
                    config=fold_cfg,
                    run_id=run_id,
                    device_str=args.device,
                    max_epochs=max_epochs,
                    max_loci=max_loci,
                )
            elif kind == "hier_train":
                result = train_hierarchical_baseline(
                    project_root=paths.project_root,
                    data_root=paths.data_root,
                    artifact_root=paths.artifact_root,
                    config=fold_cfg,
                    run_id=run_id,
                    device_str=args.device,
                    max_epochs=max_epochs,
                    max_loci=max_loci,
                )
            else:
                raise ValueError(f"unsupported arm kind: {kind!r}")
            metrics = result.metrics
            fold_rows.append({"fold": fold_idx, "run_id": run_id, "metrics": metrics, "resumed": False})

        f1_vals = [_tissue_f1(r["metrics"]) for r in fold_rows]
        nums = [v for v in f1_vals if v is not None]
        results["arms"][arm_id] = {
            "kind": kind,
            "folds": fold_rows,
            "mean_tissue_f1": float(np.mean(nums)) if nums else None,
        }
        write_json(report_dir / "per_arm" / f"{arm_id}.json", results["arms"][arm_id])

    write_json(report_dir / "summary.json", results)
    print(f"[encoder-parity] wrote {report_dir / 'summary.json'}", flush=True)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(130)
