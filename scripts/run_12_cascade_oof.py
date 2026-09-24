#!/usr/bin/env python3
"""Milestone 12 -- P2-G cascade 5x6 OOF (second arm, after N-light).

Trains the cascade encoder (P2-G scalar max/max -- the confirmed nine-pack
finalist; the 10e staged S1-S4 recipe was tried and rejected, see
milestone-7h-pretrained-mbs-rbs-campaign.md 2026-09-08) on
``hub-nine-pack-5fold-v1``. Joint multitask heads (age/tissue/sex) are
encoder supervision only. Product readout is post-hoc ``mbs_enet_nested`` /
``rbs_enet_nested`` (CPU, overlaps the next GPU fold), matching the N-light
OOF runner's pattern.

Skip-if-done on ``metrics.json``. One (fold, restart) combination per
``run_cascade_hub`` call via ``fold_indices``/``seed_override`` (added
2026-09-08 specifically for this runner -- ``run_cascade_hub`` previously
only supported a "first N folds" prefix slice, not arbitrary fold selection
or restart-independent seeds).
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from copy import deepcopy
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mbs.paths import DataPaths  # noqa: E402
from mbs.training.cascade_loop import run_cascade_hub  # noqa: E402
from mbs.training.dev_cv import load_frozen_folds  # noqa: E402
from mbs.training.loop import load_experiment_config  # noqa: E402

DEFAULT_CONFIG = ROOT / "configs" / "experiment" / "stage0_12_cascade_oof.yaml"
DEFAULT_PREFIX = "stage0-12-cascade-oof"
DEFAULT_REPORT = ROOT / "reports" / "inspection" / "stage0_12_cascade_oof"
RESTART_SEEDS = (42, 43, 44, 45, 46, 47)


def _run_id(prefix: str, fold_i: int, restart: int) -> str:
    return f"{prefix}-f{fold_i}-r{restart}"


def _has_nested(metrics: dict[str, Any]) -> bool:
    ev = metrics.get("evaluations") or {}
    return isinstance(ev.get("mbs_enet_nested"), dict)


def _enqueue_nested(
    *,
    config_path: Path,
    run_id: str,
    executor: ThreadPoolExecutor,
    pending: list,
) -> None:
    cmd = [
        sys.executable,
        str(ROOT / "scripts" / "eval_mbs_enet_from_scores.py"),
        "--run-id",
        run_id,
        "--config",
        str(config_path),
        "--which",
        "both",
        "--nested",
    ]
    print(f"[12-cascade] queue nested enet {run_id}", flush=True)
    pending.append(
        (run_id, executor.submit(subprocess.run, cmd, check=False, cwd=str(ROOT)))
    )


def _drain(pending: list, *, block: bool = False) -> None:
    still = []
    for run_id, fut in pending:
        if block or fut.done():
            result = fut.result()
            code = getattr(result, "returncode", 0)
            print(f"[12-cascade] nested enet {run_id} exit={code}", flush=True)
        else:
            still.append((run_id, fut))
    pending[:] = still


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--run-prefix", default=DEFAULT_PREFIX)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--folds", default="all", help="all or comma-separated indices")
    parser.add_argument("--restarts", default="all", help="all or comma-separated 0..5")
    parser.add_argument("--skip-enet", action="store_true")
    args = parser.parse_args()

    paths = DataPaths.from_environment()
    cfg = load_experiment_config(args.config)
    split_id = str(cfg.get("split_id", "hub-nine-pack-5fold-v1"))
    fold_pack = load_frozen_folds(paths.artifact_root / "splits" / split_id / "folds.json")
    n_folds = len(fold_pack["folds"])
    fold_idx = (
        list(range(n_folds))
        if args.folds == "all"
        else [int(x) for x in str(args.folds).split(",") if x.strip() != ""]
    )
    restart_idx = (
        list(range(len(RESTART_SEEDS)))
        if args.restarts == "all"
        else [int(x) for x in str(args.restarts).split(",") if x.strip() != ""]
    )
    report_dir = args.report_dir if args.report_dir.is_absolute() else ROOT / args.report_dir
    report_dir.mkdir(parents=True, exist_ok=True)

    pending: list = []
    executor = ThreadPoolExecutor(max_workers=1)
    jobs = [(f, r) for f in fold_idx for r in restart_idx]
    print(
        f"[12-cascade] split={split_id} n_folds={n_folds} jobs={len(jobs)} "
        f"device={args.device} prefix={args.run_prefix}",
        flush=True,
    )
    try:
        for fold_i, restart in jobs:
            if fold_i < 0 or fold_i >= n_folds:
                raise SystemExit(f"fold {fold_i} out of range 0..{n_folds - 1}")
            if restart < 0 or restart >= len(RESTART_SEEDS):
                raise SystemExit(f"restart {restart} out of range 0..{len(RESTART_SEEDS) - 1}")
            run_id = _run_id(args.run_prefix, fold_i, restart)
            run_root = paths.artifact_root / "runs" / run_id
            # run_cascade_hub writes per-fold metrics under runs/<run_id>/fold_<i>/.
            metrics_path = run_root / f"fold_{fold_i}" / "metrics.json"
            seed = int(RESTART_SEEDS[restart])
            if metrics_path.is_file():
                metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
                print(f"[12-cascade] skip-if-done train {run_id}", flush=True)
                if not args.skip_enet and not _has_nested(metrics):
                    _enqueue_nested(
                        config_path=args.config, run_id=run_id, executor=executor, pending=pending
                    )
                _drain(pending)
                continue

            fold_cfg = deepcopy(cfg)
            print(
                f"[12-cascade] train {run_id} fold={fold_i} restart={restart} seed={seed}",
                flush=True,
            )
            run_cascade_hub(
                project_root=paths.project_root,
                data_root=paths.data_root,
                artifact_root=paths.artifact_root,
                config=fold_cfg,
                run_id=run_id,
                device_str=args.device,
                fold_indices=[fold_i],
                seed_override=seed,
                report_dir=report_dir / f"_staging_{run_id.replace('-', '_')}",
                skip_if_done=True,
            )
            if not args.skip_enet:
                _enqueue_nested(
                    config_path=args.config, run_id=run_id, executor=executor, pending=pending
                )
            _drain(pending)
        _drain(pending, block=True)
    finally:
        executor.shutdown(wait=True)

    rows = []
    for fold_i, restart in jobs:
        run_id = _run_id(args.run_prefix, fold_i, restart)
        metrics_path = paths.artifact_root / "runs" / run_id / f"fold_{fold_i}" / "metrics.json"
        if not metrics_path.is_file():
            continue
        metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
        ev = metrics.get("evaluations") or {}
        nested = (ev.get("mbs_enet_nested") or {}).get("metrics") or {}
        e2e = (ev.get("mbs_e2e") or {}).get("metrics") or {}
        rows.append(
            {
                "run_id": run_id,
                "fold": fold_i,
                "restart": restart,
                "nested_tissue_f1": (nested.get("tissue") or {}).get("macro_f1"),
                "nested_age_mae": (nested.get("age") or {}).get("mae"),
                "nested_sex_auroc": (nested.get("sex") or {}).get("auroc"),
                "e2e_tissue_f1": (e2e.get("tissue") or {}).get("macro_f1"),
                "e2e_age_mae": (e2e.get("age") or {}).get("mae"),
                "e2e_sex_auroc": (e2e.get("sex") or {}).get("auroc"),
            }
        )
    summary = {
        "split_id": split_id,
        "n_folds": n_folds,
        "n_restarts": len(RESTART_SEEDS),
        "run_prefix": args.run_prefix,
        "primary_evaluation": "mbs_enet_nested",
        "rows": rows,
    }
    (report_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"[12-cascade] wrote {report_dir / 'summary.json'} n_done={len(rows)}", flush=True)


if __name__ == "__main__":
    main()
