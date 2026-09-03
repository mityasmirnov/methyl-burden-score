#!/usr/bin/env python3
"""Re-evaluate N-light Stage A mbs_e2e from saved checkpoints (no encoder retrain)."""

from __future__ import annotations

import argparse
import json
from copy import deepcopy
from pathlib import Path

from mbs.paths import DataPaths
from mbs.training.dev_cv import inject_fold_into_config, load_frozen_folds
from mbs.training.flat_stage_a_eval import complete_flat_stage_a_cpu_probes
from mbs.training.loop import load_experiment_config, train_flat_baseline


def _needs_reeval(run_root: Path) -> bool:
    manifest_path = run_root / "score_manifest.json"
    metrics_path = run_root / "metrics.json"
    if not metrics_path.is_file():
        return False
    manifest: dict = {}
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    contract = str(manifest.get("orientation_contract_version", "1"))
    if contract == "2":
        metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
        e2e = (metrics.get("evaluations") or {}).get("mbs_e2e") or {}
        if e2e.get("eval_split") == "test" and contract == "2":
            return False
    return True


def reeval_run(
    *,
    paths: DataPaths,
    config: dict,
    fold: dict,
    run_id: str,
    device: str,
    refresh_probes: bool,
) -> None:
    run_root = paths.artifact_root / "runs" / run_id
    fold_cfg = inject_fold_into_config(deepcopy(config), fold, seed=42)
    fold_cfg.setdefault("training", {})
    fold_cfg["training"]["reeval_only"] = True
    fold_cfg["training"]["max_epochs"] = 0
    fold_cfg["training"]["stage_a_defer_cpu_probes"] = False
    train_flat_baseline(
        project_root=paths.project_root,
        data_root=paths.data_root,
        artifact_root=paths.artifact_root,
        config=fold_cfg,
        run_id=run_id,
        device_str=device,
        max_epochs=0,
        max_loci=int((fold_cfg.get("cv_budget") or {}).get("max_loci", 65536)),
    )
    if refresh_probes:
        complete_flat_stage_a_cpu_probes(run_root)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--arm", action="append", default=["N-light-gene-mean", "N-light-gene-max"])
    parser.add_argument("--fold", type=int, default=None, help="Single fold index (default: all)")
    parser.add_argument("--refresh-probes", action="store_true")
    parser.add_argument("--force", action="store_true", help="Re-eval even if contract v2 present")
    args = parser.parse_args()

    paths = DataPaths.from_environment()
    config_path = args.config if args.config.is_absolute() else paths.project_root / args.config
    cfg = load_experiment_config(config_path)
    split_id = str(cfg.get("split_id", "hub-ats-7e-3fold-v1"))
    fold_pack = load_frozen_folds(paths.artifact_root / "splits" / split_id / "folds.json")
    arms_cfg = {str(a["id"]): a for a in (cfg.get("arms") or []) if isinstance(a, dict)}

    for arm_id in args.arm:
        arm = arms_cfg.get(arm_id)
        if arm is None:
            prefix = arm_id
        else:
            prefix = str(arm.get("run_prefix") or arm_id)
        folds = list(enumerate(fold_pack["folds"]))
        if args.fold is not None:
            folds = [(args.fold, fold_pack["folds"][args.fold])]
        for fold_i, fold in folds:
            run_id = f"{prefix}-f{fold_i}"
            run_root = paths.artifact_root / "runs" / run_id
            if not args.force and not _needs_reeval(run_root):
                print(f"[reeval] skip {run_id} (orientation contract v2 + valid e2e)", flush=True)
                continue
            print(f"[reeval] {run_id}", flush=True)
            reeval_run(
                paths=paths,
                config=cfg,
                fold=fold,
                run_id=run_id,
                device=args.device,
                refresh_probes=args.refresh_probes,
            )


if __name__ == "__main__":
    main()
