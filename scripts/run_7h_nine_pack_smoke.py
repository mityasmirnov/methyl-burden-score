#!/usr/bin/env python3
"""7H Phase 3: nine-pack smoke + full-cohort P2-G / m-only references.

Uses frozen split ``hub-nine-pack-3fold-v1`` (34,234 samples). All HM450 —
no cross-platform claim. Cascade/flat open the virtual multi-store via
``open_betas_for_matrix``.
"""

from __future__ import annotations

import argparse
import json
import statistics
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "reports" / "inspection" / "stage0_7h_nine_pack_smoke"
P2_CFG = ROOT / "configs" / "experiment" / "stage0_7h_nine_pack_p2_g.yaml"
M_ONLY_CFG = ROOT / "configs" / "experiment" / "stage0_7h_nine_pack_m_only.yaml"
P2_RUN = "stage0-7h-nine-pack-P2-G"
M_ONLY_PREFIX = "stage0-7h-nine-pack-m-only"


def _utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def _run(cmd: list[str]) -> None:
    print(f"[nine-pack] {_utc()} RUN {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, cwd=ROOT, check=True)


def loader_self_check(*, n_cols: int = 64, n_rows: int = 8) -> dict[str, Any]:
    import numpy as np
    import pandas as pd

    from mbs.matrix.virtual_hub_store import VIRTUAL_MATRIX_ID, open_betas_for_matrix

    root = ROOT / "data" / "canonical" / "matrices" / VIRTUAL_MATRIX_ID
    betas = open_betas_for_matrix(root)
    route = pd.read_parquet(root / "route.parquet")
    rows = np.arange(min(n_rows, betas.shape[0]), dtype=np.int64)
    block = np.asarray(betas[rows, :n_cols], dtype=np.float32)
    finite = float(np.isfinite(block).mean())
    packs = (
        route.loc[route["row_index"].isin(rows.tolist()), "matrix_id"]
        .astype(str)
        .value_counts()
        .to_dict()
    )
    out = {
        "matrix_id": VIRTUAL_MATRIX_ID,
        "shape": list(betas.shape),
        "checked_rows": int(rows.size),
        "checked_cols": int(n_cols),
        "finite_fraction": finite,
        "packs_in_check": packs,
        "platform_claim": "HM450_only_no_cross_platform",
    }
    print(f"[nine-pack] loader_ok {json.dumps(out)}", flush=True)
    if finite < 0.5:
        raise SystemExit(f"loader self-check: low finite fraction {finite}")
    return out


def cascade_train(
    *,
    run_id: str,
    max_folds: int | None,
    max_epochs: int | None,
    report_dir: Path,
) -> None:
    cmd = [
        "uv",
        "run",
        "mbs",
        "train",
        "cascade",
        "--config",
        str(P2_CFG),
        "--run-id",
        run_id,
        "--device",
        "cuda",
        "--skip-if-done",
        "--report-dir",
        str(report_dir),
    ]
    if max_folds is not None:
        cmd.extend(["--max-folds", str(max_folds)])
    if max_epochs is not None:
        cmd.extend(["--max-epochs", str(max_epochs)])
    _run(cmd)


def flat_train(
    *,
    fold_filter: int | None,
    max_epochs: int | None,
    run_prefix: str | None = None,
) -> None:
    sys.path.insert(0, str(ROOT / "scripts"))
    from run_7g_gene_only_probe import train_flat_region_arm  # type: ignore
    from mbs.paths import DataPaths
    from mbs.training.loop import load_experiment_config

    paths = DataPaths.from_environment()
    cfg_path = M_ONLY_CFG
    prefix = run_prefix or M_ONLY_PREFIX
    if max_epochs is not None:
        import yaml

        cfg = load_experiment_config(cfg_path)
        cfg.setdefault("cv_budget", {})["max_epochs"] = int(max_epochs)
        cfg.setdefault("training", {})["max_epochs"] = int(max_epochs)
        tmp = REPORT / f"_tmp_m_only_epochs_{max_epochs}.yaml"
        REPORT.mkdir(parents=True, exist_ok=True)
        tmp.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")
        cfg_path = tmp
    train_flat_region_arm(
        paths=paths,
        config_path=cfg_path,
        run_prefix=prefix,
        device="cuda",
        report_dir=REPORT,
        fold_filter=fold_filter,
    )


def _e2e_metrics(metrics_path: Path) -> dict[str, float | None]:
    blob = json.loads(metrics_path.read_text(encoding="utf-8"))
    ev = (blob.get("evaluations") or {}).get("mbs_e2e") or {}
    m = ev.get("metrics") or {}
    return {
        "eval_split": ev.get("eval_split"),
        "tissue_f1": (m.get("tissue") or {}).get("macro_f1"),
        "age_mae": (m.get("age") or {}).get("mae"),
        "sex_auroc": (m.get("sex") or {}).get("auroc"),
        "best_epoch": blob.get("best_epoch")
        or (blob.get("checkpoint_selection") or {}).get("best_epoch"),
    }


def collect_arm(kind: str, *, smoke: bool = False) -> dict[str, Any]:
    folds: list[dict[str, Any]] = []
    if kind == "p2":
        run = f"{P2_RUN}-smoke" if smoke else P2_RUN
        root = ROOT / "artifacts" / "runs" / run
        for i in range(3):
            mp = root / f"fold_{i}" / "metrics.json"
            if mp.is_file():
                row = _e2e_metrics(mp)
                row["fold"] = i
                row["run_id"] = run
                folds.append(row)
    else:
        prefix = f"{M_ONLY_PREFIX}-smoke" if smoke else M_ONLY_PREFIX
        for i in range(3):
            mp = ROOT / "artifacts" / "runs" / f"{prefix}-f{i}" / "metrics.json"
            if mp.is_file():
                row = _e2e_metrics(mp)
                row["fold"] = i
                row["run_id"] = f"{prefix}-f{i}"
                folds.append(row)
        # Legacy: smoke once wrote into the full prefix (fold 0 only). Prefer
        # full-budget metrics when both exist; otherwise keep legacy for reports.
        if smoke and not folds:
            mp = ROOT / "artifacts" / "runs" / f"{M_ONLY_PREFIX}-f0" / "metrics.json"
            if mp.is_file():
                row = _e2e_metrics(mp)
                row["fold"] = 0
                row["run_id"] = f"{M_ONLY_PREFIX}-f0"
                folds.append(row)
    def _mean(key: str) -> float | None:
        vals = [float(f[key]) for f in folds if isinstance(f.get(key), (int, float))]
        return float(statistics.mean(vals)) if vals else None

    def _sd(key: str) -> float | None:
        vals = [float(f[key]) for f in folds if isinstance(f.get(key), (int, float))]
        return float(statistics.stdev(vals)) if len(vals) >= 2 else None

    return {
        "arm": "P2-G" if kind == "p2" else "N-light-gene-ablation-m-only",
        "n_folds": len(folds),
        "folds": folds,
        "tissue_f1_mean": _mean("tissue_f1"),
        "tissue_f1_sd": _sd("tissue_f1"),
        "age_mae_mean": _mean("age_mae"),
        "age_mae_sd": _sd("age_mae"),
        "sex_auroc_mean": _mean("sex_auroc"),
        "sex_auroc_sd": _sd("sex_auroc"),
    }


def write_report(loader: dict[str, Any], *, smoke: bool) -> None:
    REPORT.mkdir(parents=True, exist_ok=True)
    p2 = collect_arm("p2", smoke=smoke)
    mo = collect_arm("m_only", smoke=smoke)
    summary = {
        "generated_at": _utc(),
        "split_id": "hub-nine-pack-3fold-v1",
        "matrix_id": "matrix-hub-nine-pack-virtual-v1",
        "n_samples": 34234,
        "platform": "HM450",
        "cross_platform_claim": False,
        "phase": "smoke" if smoke else "full_reference",
        "loader_self_check": loader,
        "arms": {"P2-G": p2, "m_only": mo},
    }
    (REPORT / "summary.json").write_text(
        json.dumps(summary, indent=2, default=str) + "\n", encoding="utf-8"
    )
    lines = [
        "# Nine-pack P2-G / m-only reference snapshot",
        "",
        f"Generated: `{summary['generated_at']}`",
        "",
        "Auto-exported pair only. Curated campaign board with pooling grid,",
        "ablations, and interpretations: [`analysis.md`](analysis.md).",
        "",
        f"- Matrix: `{summary['matrix_id']}` (virtual multi-store)",
        f"- Split: `{summary['split_id']}` (**{summary['n_samples']}** samples)",
        f"- Platform: **{summary['platform']} only** — no cross-platform claim",
        f"- Phase: `{summary['phase']}`",
        "",
        "## Loader self-check",
        "",
        f"`{json.dumps(loader)}`",
        "",
        "## Results (`mbs_e2e`, outer test)",
        "",
        "| Arm | folds | Tissue F1 | Age MAE | Sex AUROC |",
        "|-----|------:|----------:|--------:|----------:|",
    ]
    for label, arm in (("P2-G cascade", p2), ("one-hop m-only", mo)):
        t = arm.get("tissue_f1_mean")
        a = arm.get("age_mae_mean")
        s = arm.get("sex_auroc_mean")
        td = arm.get("tissue_f1_sd")
        ad = arm.get("age_mae_sd")
        sd = arm.get("sex_auroc_sd")

        def fmt(m: float | None, d: float | None) -> str:
            if m is None:
                return "—"
            if d is None:
                return f"{m:.3f}"
            return f"{m:.3f} (±{d:.3f})"

        lines.append(
            f"| {label} | {arm['n_folds']} | {fmt(t, td)} | {fmt(a, ad)} | {fmt(s, sd)} |"
        )
    lines.extend(
        [
            "",
            "### Fold detail",
            "",
            "```json",
            json.dumps({"P2-G": p2.get("folds"), "m_only": mo.get("folds")}, indent=2),
            "```",
            "",
            "## Notes",
            "",
            "- Full cohort (no subsample).",
            "- 7H remains open after this smoke; trait expansion still gated on ≥1k-per-arm census.",
            "- Stage B / Milestone 7 OOF stay blocked.",
            "",
        ]
    )
    # Keep analysis.md as the curated campaign board (pooling grid, ablations,
    # interpretations). This auto snapshot is only the P2-G / m-only pair.
    snapshot = REPORT / "p2_m_only_reference.md"
    snapshot.write_text("\n".join(lines), encoding="utf-8")
    print(f"[nine-pack] wrote {snapshot} (analysis.md is curated; not overwritten)", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--phase",
        choices=("check", "smoke", "full", "report"),
        default="full",
        help="check=loader only; smoke=fold0 short; full=3-fold refs; report=rewrite metrics",
    )
    parser.add_argument("--skip-p2", action="store_true")
    parser.add_argument("--skip-m-only", action="store_true")
    args = parser.parse_args()

    REPORT.mkdir(parents=True, exist_ok=True)
    loader = loader_self_check()
    (REPORT / "loader_self_check.json").write_text(
        json.dumps(loader, indent=2) + "\n", encoding="utf-8"
    )
    if args.phase == "check":
        return

    if args.phase == "smoke":
        if not args.skip_p2:
            cascade_train(
                run_id=f"{P2_RUN}-smoke",
                max_folds=1,
                max_epochs=3,
                report_dir=REPORT / "_staging_p2_smoke",
            )
        if not args.skip_m_only:
            # Distinct run prefix so smoke never poisons full skip-if-done.
            flat_train(
                fold_filter=0,
                max_epochs=3,
                run_prefix=f"{M_ONLY_PREFIX}-smoke",
            )
        write_report(loader, smoke=True)
        return

    if args.phase == "full":
        if not args.skip_p2:
            cascade_train(
                run_id=P2_RUN,
                max_folds=None,
                max_epochs=None,
                report_dir=REPORT / "_staging_p2_g",
            )
        if not args.skip_m_only:
            flat_train(fold_filter=None, max_epochs=None)
        write_report(loader, smoke=False)
        return

    write_report(loader, smoke=False)


if __name__ == "__main__":
    main()
