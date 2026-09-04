#!/usr/bin/env python3
"""Compute Stage A ablation representation diagnostics from saved MBS + checkpoints."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from mbs.annotation.manifest import write_json
from mbs.inspection.repr_diagnostics import (
    cache_sample_mean_m_gene_panel,
    compute_mbs_repr_diagnostics,
    head_weight_l2_from_checkpoint,
)
from mbs.paths import DataPaths

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REPORT = ROOT / "reports/inspection/stage0_7g_gene_only_probe"

# short label -> per_arm arm_id base (seeds: base + base-s2)
ABLATION_ARMS = [
    ("A0", "N-light-gene-ablation-m-only"),
    ("A1", "N-light-gene-ablation-m-role"),
    ("A2", "N-light-gene-ablation-m-context"),
    ("A3", "N-light-gene-ablation-m-role-context"),
    ("A4/A7", "N-light-gene-ablation-full"),
    ("N0", "N-light-gene-ablation-n0-obs-only"),
    ("N1", "N-light-gene-ablation-n1-anno-only"),
    ("N2", "N-light-gene-ablation-n2-reg-permuted"),
    ("N3", "N-light-gene-ablation-n3-reg-zero"),
]


def _seed_ids(base: str) -> list[str]:
    return [base, f"{base}-s2"]


def _run_from_fold(fold: dict[str, Any]) -> Path | None:
    score_dir = fold.get("score_dir")
    if not score_dir:
        return None
    return Path(str(score_dir)).parent


def _ckpt_for_run(run_dir: Path) -> Path | None:
    name = run_dir.name
    cand = [
        run_dir.parent.parent / "checkpoints" / name / "best.pt",
        ROOT / "artifacts" / "checkpoints" / name / "best.pt",
        run_dir / "best.pt",
    ]
    for p in cand:
        if p.is_file():
            return p
    return None


def _diagnostics_for_fold(
    fold: dict[str, Any],
    *,
    sample_mean_m: np.ndarray | None,
) -> dict[str, Any] | None:
    run_dir = _run_from_fold(fold)
    if run_dir is None:
        return None
    mbs_path = run_dir / "scores" / "mbs.npy"
    present_path = run_dir / "scores" / "mbs_present.npy"
    if not mbs_path.is_file():
        return None
    mbs = np.load(mbs_path)
    present = np.load(present_path) if present_path.is_file() else None
    ckpt = _ckpt_for_run(run_dir)
    head_l2 = head_weight_l2_from_checkpoint(ckpt) if ckpt else None
    sel = fold.get("checkpoint_selection") or {}
    best_epoch = sel.get("best_epoch") or fold.get("best_epoch")
    best_val = None
    hist = sel.get("val_history") or []
    if best_epoch is not None and hist:
        for row in hist:
            if int(row.get("epoch", -1)) == int(best_epoch):
                best_val = row.get("tissue_macro_f1")
                break
    diag = compute_mbs_repr_diagnostics(
        mbs,
        present,
        sample_mean_m=sample_mean_m,
        head_weight_l2=head_l2,
        best_epoch=int(best_epoch) if best_epoch is not None else None,
        best_val_tissue_f1=float(best_val) if best_val is not None else None,
    )
    diag["run_dir"] = str(run_dir)
    diag["checkpoint"] = str(ckpt) if ckpt else None
    return diag


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--skip-mean-m", action="store_true")
    args = parser.parse_args()
    report_dir: Path = args.report_dir
    per_arm = report_dir / "per_arm"
    paths = DataPaths.from_environment()

    sample_mean_m = None
    if not args.skip_mean_m:
        panel = json.loads((report_dir / "gene_panel_manifest.json").read_text(encoding="utf-8"))
        cols = np.asarray(panel["gene_col_indices"], dtype=np.int64)
        betas = (
            paths.data_root
            / "canonical"
            / "matrices"
            / str(panel["matrix_id"])
            / "betas.zarr"
        )
        cache_path = report_dir / "sample_mean_m_gene_panel.npy"
        print(f"[repr] caching sample mean-M → {cache_path}", flush=True)
        sample_mean_m = cache_sample_mean_m_gene_panel(
            betas_zarr=betas,
            gene_col_indices=cols,
            out_path=cache_path,
        )

    payload: dict[str, Any] = {"arms": {}, "by_label": {}}
    for label, base in ABLATION_ARMS:
        seed_diags: list[dict[str, Any]] = []
        for arm_id in _seed_ids(base):
            path = per_arm / f"{arm_id}.json"
            if not path.is_file():
                continue
            arm = json.loads(path.read_text(encoding="utf-8"))
            folds = arm.get("folds") or []
            fold_diags = []
            for fold in folds:
                d = _diagnostics_for_fold(fold, sample_mean_m=sample_mean_m)
                if d is None:
                    continue
                fold_diags.append(d)
                ev = fold.setdefault("evaluations", {})
                e2e = ev.setdefault("mbs_e2e", {})
                e2e["repr_diagnostics"] = {
                    k: d[k]
                    for k in (
                        "gene_score_sd",
                        "mean_per_gene_sd",
                        "saturation_fraction",
                        "constant_score_fraction",
                        "corr_mean_m",
                        "head_weight_l2",
                        "best_epoch",
                        "best_val_tissue_f1",
                        "score_mean",
                        "score_min",
                        "score_max",
                    )
                    if k in d
                }
            if not fold_diags:
                continue
            path.write_text(json.dumps(arm) + "\n", encoding="utf-8")
            mean = {
                k: float(np.mean([fd[k] for fd in fold_diags if k in fd and np.isfinite(fd[k])]))
                for k in (
                    "gene_score_sd",
                    "mean_per_gene_sd",
                    "saturation_fraction",
                    "constant_score_fraction",
                    "corr_mean_m",
                    "head_weight_l2",
                    "best_val_tissue_f1",
                )
            }
            mean["best_epoch"] = float(
                np.mean([fd["best_epoch"] for fd in fold_diags if "best_epoch" in fd])
            ) if any("best_epoch" in fd for fd in fold_diags) else None
            payload["arms"][arm_id] = {"folds": fold_diags, "mean": mean}
            seed_diags.append(mean)
        if seed_diags:
            payload["by_label"][label] = {
                "mean": {
                    k: float(np.mean([s[k] for s in seed_diags if s.get(k) is not None and np.isfinite(s[k])]))
                    for k in (
                        "gene_score_sd",
                        "mean_per_gene_sd",
                        "saturation_fraction",
                        "constant_score_fraction",
                        "corr_mean_m",
                        "head_weight_l2",
                        "best_val_tissue_f1",
                    )
                }
            }
            print(
                f"[repr] {label}: sd={payload['by_label'][label]['mean'].get('gene_score_sd'):.4f} "
                f"sat={payload['by_label'][label]['mean'].get('saturation_fraction'):.4f} "
                f"const={payload['by_label'][label]['mean'].get('constant_score_fraction'):.4f} "
                f"corr_m={payload['by_label'][label]['mean'].get('corr_mean_m'):.4f}",
                flush=True,
            )

    out = report_dir / "repr_diagnostics.json"
    write_json(out, payload)
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
