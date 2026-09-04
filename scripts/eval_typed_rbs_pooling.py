#!/usr/bin/env python3
"""Evaluate typed-RBS pooling arms (R0..R5 + role-shuffle control) on CPU.

Reads the per-fold ``all_gene_rbs`` region matrix saved by a gene-only probe run
and scores each pooling arm with the transparent linear multitask heads. Writes
a report under ``reports/inspection/stage0_7g_gene_only_probe/typed_rbs_pooling/``.

Example::

    source scripts/activate_data_environment.sh
    python scripts/eval_typed_rbs_pooling.py \
        --run-id stage0-7g-gene-probe-vector-mean-max
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import zarr

from mbs.annotation.manifest import utc_now_iso, write_json
from mbs.paths import DataPaths
from mbs.training.dev_cv import load_frozen_folds
from mbs.training.loop import load_experiment_config
from mbs.training.phenotypes import load_multitask_phenotypes
from mbs.training.typed_rbs_pooling import (
    arm_builders,
    evaluate_arm,
    features_typed,
    shuffle_region_types,
    typed_pool_promotion_gate,
)

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "configs/experiment/stage0_7g_gene_only_probe.yaml"
DEFAULT_RUN = "stage0-7g-gene-probe-vector-mean-max"
REPORT_DIR = ROOT / "reports/inspection/stage0_7g_gene_only_probe/typed_rbs_pooling"


def _fold_index(fold: dict, sample_ids: list[str]) -> tuple[np.ndarray, np.ndarray]:
    train_set = {str(s) for s in fold["train_sample_ids"]}
    test_ids = fold.get("external_test_sample_ids") or fold.get("validation_sample_ids") or []
    test_set = {str(s) for s in test_ids}
    train_idx = np.asarray(
        [i for i, sid in enumerate(sample_ids) if sid in train_set], dtype=np.int64
    )
    test_idx = np.asarray(
        [i for i, sid in enumerate(sample_ids) if sid in test_set], dtype=np.int64
    )
    if train_idx.size == 0 or test_idx.size == 0:
        raise ValueError("empty train or test index")
    return train_idx, test_idx


def _phenotype_arrays(sample_ids: list[str], ph_by_id: dict) -> dict[str, np.ndarray]:
    return {
        "age": np.asarray([float(ph_by_id[s].age or 0.0) for s in sample_ids], dtype=np.float64),
        "age_mask": np.asarray([bool(ph_by_id[s].age_mask) for s in sample_ids], dtype=bool),
        "tissue": np.asarray(
            [int(ph_by_id[s].class_index) if ph_by_id[s].tissue_mask else 0 for s in sample_ids],
            dtype=np.int64,
        ),
        "tissue_mask": np.asarray([bool(ph_by_id[s].tissue_mask) for s in sample_ids], dtype=bool),
        "sex": np.asarray(
            [
                int(ph_by_id[s].sex_class_index or 0) if ph_by_id[s].sex_mask else 0
                for s in sample_ids
            ],
            dtype=np.int64,
        ),
        "sex_mask": np.asarray([bool(ph_by_id[s].sex_mask) for s in sample_ids], dtype=bool),
        "study_ids": np.asarray(
            [str(ph_by_id[s].study_id or "NA") for s in sample_ids], dtype=object
        ),
    }


def _slice(arrays: dict[str, np.ndarray], idx: np.ndarray) -> dict[str, np.ndarray]:
    return {k: v[idx] for k, v in arrays.items()}


def _scalars(result: dict[str, Any]) -> dict[str, float | None]:
    m = result.get("metrics") or {}
    return {
        "age_mae": (m.get("age") or {}).get("mae"),
        "age_r2": (m.get("age") or {}).get("r2"),
        "tissue_f1": (m.get("tissue") or {}).get("macro_f1"),
        "sex_auroc": (m.get("sex") or {}).get("auroc"),
        "n_features": result.get("n_features"),
    }


def _mean_scalars(per_fold: list[dict[str, float | None]]) -> dict[str, float | None]:
    keys = ("age_mae", "age_r2", "tissue_f1", "sex_auroc", "n_features")
    out: dict[str, float | None] = {}
    for k in keys:
        vals = [float(f[k]) for f in per_fold if f.get(k) is not None]
        out[k] = float(np.mean(vals)) if vals else None
    return out


def evaluate_run(
    *,
    run_root: Path,
    folds: list[dict],
    ph_by_id: dict[str, Any],
    class_names: list[str],
    include_r5: bool,
    control_seeds: list[int],
) -> dict[str, Any]:
    builders = arm_builders(include_r5=include_r5)
    per_arm_folds: dict[str, list[dict[str, float | None]]] = {a: [] for a in builders}
    control_folds: list[dict[str, float | None]] = []
    control_stats: list[dict[str, float]] = []

    for i, fold in enumerate(folds):
        score_dir = run_root / f"fold_{i}" / "scores"
        sample_index = pd.read_parquet(score_dir / "sample_index.parquet")
        sample_ids = sample_index.sort_values("row_index")["sample_id"].astype(str).tolist()
        region_df = pd.read_parquet(score_dir / "all_gene_region_index.parquet")
        rbs = np.asarray(
            zarr.open_array(str(score_dir / "all_gene_rbs.zarr"), mode="r"), dtype=np.float32
        )
        present = np.asarray(
            zarr.open_array(str(score_dir / "all_gene_rbs_present.zarr"), mode="r"), dtype=np.uint8
        )
        train_idx, test_idx = _fold_index(fold, sample_ids)
        arrays = _phenotype_arrays(sample_ids, ph_by_id)
        ph_tr, ph_te = _slice(arrays, train_idx), _slice(arrays, test_idx)

        for arm, builder in builders.items():
            feats, _names = builder(rbs, present, region_df)
            result = evaluate_arm(
                name=arm,
                x_train=feats[train_idx],
                x_test=feats[test_idx],
                pheno_train=ph_tr,
                pheno_test=ph_te,
                class_names=class_names,
            )
            per_arm_folds[arm].append(_scalars(result))
            print(f"[{arm}] fold{i} {_scalars(result)}", flush=True)

        # Control: role-shuffle within gene, then R1.
        for seed in control_seeds:
            shuffled_df, stats = shuffle_region_types(region_df, seed=seed)
            feats, _names = features_typed(rbs, present, shuffled_df, stats=("max",))
            result = evaluate_arm(
                name=f"control_seed{seed}",
                x_train=feats[train_idx],
                x_test=feats[test_idx],
                pheno_train=ph_tr,
                pheno_test=ph_te,
                class_names=class_names,
            )
            sc = _scalars(result)
            control_folds.append(sc)
            control_stats.append(stats)
            print(f"[control seed{seed}] fold{i} {sc} {stats}", flush=True)

    arms = {a: _mean_scalars(f) for a, f in per_arm_folds.items()}
    control = _mean_scalars(control_folds) if control_folds else {}
    return {
        "arms": arms,
        "arms_per_fold": per_arm_folds,
        "control": control,
        "control_stats": control_stats,
    }


def _write_report(payload: dict[str, Any], gate: dict[str, Any]) -> None:
    write_json(REPORT_DIR / "summary.json", payload)
    arms = payload["arms"]
    lines = [
        "# Typed-RBS pooling (Milestone 7G)",
        "",
        f"Run: `{payload['run_id']}` · folds: {payload['n_folds']} · "
        f"generated {payload['generated_at']}",
        "",
        "## Arms (fold-mean holdout)",
        "",
        "| arm | age MAE | age R² | tissue F1 | sex AUROC | n_feat |",
        "|-----|--------:|-------:|----------:|----------:|-------:|",
    ]

    def fmt(v: float | None, nd: int = 3) -> str:
        return f"{v:.{nd}f}" if isinstance(v, (int, float)) else "—"

    for arm in sorted(arms):
        s = arms[arm]
        lines.append(
            f"| {arm} | {fmt(s['age_mae'])} | {fmt(s['age_r2'])} | {fmt(s['tissue_f1'])} "
            f"| {fmt(s['sex_auroc'])} | {int(s['n_features']) if s['n_features'] else '—'} |"
        )
    ctrl = payload.get("control") or {}
    if ctrl:
        lines.append(
            f"| control(R1,shuffled) | {fmt(ctrl.get('age_mae'))} | {fmt(ctrl.get('age_r2'))} "
            f"| {fmt(ctrl.get('tissue_f1'))} | {fmt(ctrl.get('sex_auroc'))} | — |"
        )
    cstats = payload.get("control_stats") or []
    if cstats:
        fga = float(np.mean([c["frac_genes_altered"] for c in cstats]))
        fcc = float(np.mean([c["frac_columns_changed"] for c in cstats]))
        lines += [
            "",
            "## Role-shuffle control",
            "",
            f"- permutations: {len(cstats)} (seeds × folds)",
            f"- frac genes altered: {fga:.3f}",
            f"- frac RBS columns whose role changed: {fcc:.3f}",
        ]
    lines += [
        "",
        "## Promotion gate (R1 typed vs R0 untyped)",
        "",
        f"- **promote: {gate['promote']}**",
        f"- age MAE gain: {fmt(gate['age_mae_gain'])} (≥1.0 to pass alone)",
        f"- age R² gain: {fmt(gate['age_r2_gain'])} (≥0.05 to pass alone)",
        f"- tissue F1 loss: {fmt(gate['tissue_f1_loss'])} (≤0.03)",
        f"- sex AUROC loss: {fmt(gate['sex_auroc_loss'])} (≤0.03 if present)",
        f"- shuffle worse than typed: {gate['shuffle_worse_than_typed']}",
        "",
    ]
    (REPORT_DIR / "analysis.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", default=DEFAULT_RUN)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--include-r5", action="store_true", help="Add slow no-pool ceiling arm")
    parser.add_argument(
        "--control-seeds",
        type=int,
        nargs="+",
        default=[0, 1, 2],
        help="Fixed role-shuffle permutation seeds",
    )
    args = parser.parse_args()

    paths = DataPaths.from_environment()
    cfg = load_experiment_config(args.config)
    split_id = str(cfg.get("split_id", "hub-ats-7e-3fold-v1"))
    fold_pack = load_frozen_folds(paths.artifact_root / "splits" / split_id / "folds.json")
    folds = fold_pack["folds"]
    pheno_rel = Path(
        str(
            cfg.get("sample_phenotype_table")
            or (cfg.get("pilot") or {}).get("sample_phenotype_table")
            or "canonical/phenotypes/sample_phenotype_table_age_tissue_sex_full_v1.parquet"
        )
    )
    pheno_path = pheno_rel if pheno_rel.is_absolute() else paths.data_root / pheno_rel
    phenotypes, class_names = load_multitask_phenotypes(pheno_path)
    ph_by_id = {p.sample_id: p for p in phenotypes}

    run_root = paths.artifact_root / "runs" / args.run_id
    if not (run_root / "fold_0").is_dir():
        raise FileNotFoundError(f"run not found or missing fold dirs: {run_root}")

    results = evaluate_run(
        run_root=run_root,
        folds=folds,
        ph_by_id=ph_by_id,
        class_names=class_names,
        include_r5=bool(args.include_r5),
        control_seeds=list(args.control_seeds),
    )

    r0 = results["arms"]["R0"]
    typed = results["arms"]["R1"]
    ctrl = results.get("control") or {}
    gate = typed_pool_promotion_gate(
        age_mae_r0=float(r0["age_mae"] or 0.0),
        age_mae_typed=float(typed["age_mae"] or 0.0),
        age_r2_r0=float(r0["age_r2"] or 0.0),
        age_r2_typed=float(typed["age_r2"] or 0.0),
        tissue_f1_r0=float(r0["tissue_f1"] or 0.0),
        tissue_f1_typed=float(typed["tissue_f1"] or 0.0),
        sex_auroc_r0=r0["sex_auroc"],
        sex_auroc_typed=typed["sex_auroc"],
        shuffle_age_mae=float(ctrl.get("age_mae") or typed["age_mae"] or 0.0),
    )

    payload: dict[str, Any] = {
        "run_id": args.run_id,
        "n_folds": len(folds),
        "generated_at": utc_now_iso(),
        "include_r5": bool(args.include_r5),
        "control_seeds": list(args.control_seeds),
        "promotion_gate": gate,
        **results,
    }
    _write_report(payload, gate)
    print(json.dumps({"promotion_gate": gate}, indent=2))
    print(f"[report] wrote {REPORT_DIR}")


if __name__ == "__main__":
    main()
