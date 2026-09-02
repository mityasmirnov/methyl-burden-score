#!/usr/bin/env python3
"""Milestone 7G methylation-only full evaluation report.

Assumes cascade folds under ``$MBS_ARTIFACT_ROOT/runs/<run_id>/fold_*/`` are
complete (or runs classical/transparent only). Emits ranking table without
metadata-only arms; sex AUROC + tissue OvR from neural fusion scores.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from mbs.annotation.manifest import write_json
from mbs.inspection.arm_glossary import arm_description, render_arm_glossary_section
from mbs.paths import DataPaths
from mbs.training.classical_mvalue import run_classical_mvalue
from mbs.training.dev_cv import load_frozen_folds, samples_from_phenotype_table
from mbs.training.loop import load_experiment_config
from mbs.training.transparent_hub import run_all_transparent_arms

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REPORT = ROOT / "reports/inspection/stage0_7g_methylation_eval"
DEFAULT_SPLIT = "hub-ats-7e-3fold-v1"
DEFAULT_RUN = "stage0-7g-cascade-v1"
DEFAULT_ONTOLOGY = "canonical/phenotypes/tissue_ontology_age_tissue_sex_full_v1.yaml"
MAX_LOCI = 65536

RANKING_ARMS = (
    "N-cascade-l1",
    "T-mean-gene",
    "T-mean-region",
    "T-enet",
    "C-mvalue-ridge",
    "C-mvalue-enet",
    "C-mvalue-hgb",
    "C-mvalue-sva",
)

ARM_LABELS = {arm: arm_description(arm) for arm in RANKING_ARMS}


def _mean_std(xs: list[float | None]) -> tuple[float | None, float | None]:
    vals = [float(x) for x in xs if x is not None and not (isinstance(x, float) and math.isnan(x))]
    if not vals:
        return None, None
    arr = np.asarray(vals, dtype=np.float64)
    return float(arr.mean()), float(arr.std(ddof=1)) if arr.size > 1 else 0.0


def _style_axes(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linestyle=":", alpha=0.4)


def _metrics_from_fold_blob(blob: dict[str, Any]) -> dict[str, Any]:
    """Normalize cascade / transparent / classical fold metrics."""
    m = blob.get("metrics") if isinstance(blob.get("metrics"), dict) else blob
    age = m.get("age") if isinstance(m.get("age"), dict) else None
    tissue = m.get("tissue") if isinstance(m.get("tissue"), dict) else None
    sex = m.get("sex") if isinstance(m.get("sex"), dict) else None
    excluded = (tissue or {}).get("excluded_zero_shot_test_counts") or {}
    return {
        "tissue_macro_f1": (tissue or {}).get("macro_f1"),
        "tissue_balanced_accuracy": (tissue or {}).get("balanced_accuracy"),
        "n_classes_scored": (tissue or {}).get("n_classes_scored"),
        "n_excluded_zero_shot_samples": sum(int(v) for v in excluded.values()),
        "age_mae": (age or {}).get("mae"),
        "age_r2": (age or {}).get("r2"),
        "sex_auroc": (sex or {}).get("auroc"),
        "tissue_roc": m.get("tissue_roc"),
        "sex": sex,
    }


def cascade_rows(cascade_summary: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for fold_i, fold in enumerate(cascade_summary.get("folds") or []):
        m = _metrics_from_fold_blob(fold)
        rows.append(
            {
                "arm": "N-cascade-l1",
                "family": "neural",
                "fold": fold_i,
                **m,
            }
        )
    return rows


def classical_rows(classical: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for fold in classical.get("folds") or []:
        for arm, blob in (fold.get("arms") or {}).items():
            m = _metrics_from_fold_blob(blob)
            rows.append({"arm": arm, "family": "classical", "fold": fold["fold"], **m})
    return rows


def transparent_rows(transparent: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for fold in transparent.get("folds") or []:
        for arm, blob in (fold.get("arms") or {}).items():
            m = _metrics_from_fold_blob(blob)
            rows.append({"arm": arm, "family": "transparent", "fold": fold["fold"], **m})
    return rows


def aggregate_arms(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_arm: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for r in rows:
        if r["arm"] not in RANKING_ARMS:
            continue
        by_arm[str(r["arm"])].append(r)
    out = []
    for arm in RANKING_ARMS:
        items = by_arm.get(arm) or []
        if not items:
            continue
        f1, f1s = _mean_std([r.get("tissue_macro_f1") for r in items])
        bacc, baccs = _mean_std([r.get("tissue_balanced_accuracy") for r in items])
        mae, maes = _mean_std([r.get("age_mae") for r in items])
        r2, r2s = _mean_std([r.get("age_r2") for r in items])
        sex_a, sex_s = _mean_std([r.get("sex_auroc") for r in items])
        n_scored, _ = _mean_std([r.get("n_classes_scored") for r in items])
        n_excluded_total = sum(int(r.get("n_excluded_zero_shot_samples") or 0) for r in items)
        out.append(
            {
                "arm": arm,
                "label": ARM_LABELS.get(arm, arm),
                "family": items[0].get("family"),
                "n_cells": len(items),
                "tissue_macro_f1": f1,
                "tissue_macro_f1_sd": f1s,
                "tissue_balanced_accuracy": bacc,
                "tissue_balanced_accuracy_sd": baccs,
                "tissue_n_classes_scored_mean": n_scored,
                "tissue_n_excluded_zero_shot_samples": n_excluded_total,
                "age_mae_years": mae,
                "age_mae_years_sd": maes,
                "age_r2": r2,
                "age_r2_sd": r2s,
                "sex_auroc": sex_a,
                "sex_auroc_sd": sex_s,
            }
        )
    return out


def pick_winner(means: list[dict[str, Any]]) -> dict[str, Any] | None:
    """Max tissue macro-F1, then min age MAE among methylation-input arms."""
    ranked = [
        r
        for r in means
        if r.get("tissue_macro_f1") is not None and r["arm"] != "C-metadata"
    ]
    if not ranked:
        return None
    ranked.sort(
        key=lambda r: (
            -float(r["tissue_macro_f1"]),
            float(r["age_mae_years"]) if r.get("age_mae_years") is not None else 1e9,
        )
    )
    best = dict(ranked[0])
    best["reason"] = (
        "max mean tissue macro-F1, then min mean age MAE among methylation-input methods"
    )
    return best


def plot_roc(
    curves: list[dict[str, Any]] | None,
    sex: dict[str, Any] | None,
    tissue_path: Path,
    sex_path: Path,
    *,
    title_prefix: str,
) -> None:
    if curves:
        fig, ax = plt.subplots(figsize=(6.4, 5.2))
        ax.plot([0, 1], [0, 1], linestyle="--", color="#888", linewidth=1, label="Chance")
        for curve in curves:
            ax.plot(
                curve["fpr"],
                curve["tpr"],
                label=f"{curve['label']} (AUROC {curve['auroc']:.2f})",
            )
        ax.set_xlabel("False positive rate")
        ax.set_ylabel("True positive rate")
        ax.set_title(f"{title_prefix} — tissue one-vs-rest ROC")
        ax.legend(fontsize=8, loc="lower right")
        _style_axes(ax)
        fig.tight_layout()
        fig.savefig(tissue_path, dpi=140)
        plt.close(fig)
    if isinstance(sex, dict) and sex.get("fpr") and sex.get("tpr"):
        fig, ax = plt.subplots(figsize=(5.4, 5.0))
        ax.plot([0, 1], [0, 1], linestyle="--", color="#888", linewidth=1, label="Chance")
        auroc = sex.get("auroc")
        label = f"Sex (AUROC {auroc:.2f})" if auroc is not None else "Sex"
        ax.plot(sex["fpr"], sex["tpr"], label=label)
        ax.set_xlabel("False positive rate")
        ax.set_ylabel("True positive rate")
        ax.set_title(f"{title_prefix} — sex ROC")
        ax.legend(loc="lower right")
        _style_axes(ax)
        fig.tight_layout()
        fig.savefig(sex_path, dpi=140)
        plt.close(fig)


def plot_ranking_bars(means: list[dict[str, Any]], path: Path) -> None:
    names, f1s, maes, sexes = [], [], [], []
    for row in means:
        names.append(ARM_LABELS.get(row["arm"], row["arm"]))
        f1s.append(row["tissue_macro_f1"] or float("nan"))
        mae = row["age_mae_years"]
        maes.append(mae if mae is not None and 0 < mae < 100 else float("nan"))
        sexes.append(row["sex_auroc"] if row["sex_auroc"] is not None else float("nan"))
    fig, axes = plt.subplots(1, 3, figsize=(14, 5.2))
    y = np.arange(len(names))
    axes[0].barh(y, f1s, color="#3b6d9a")
    axes[0].set_yticks(y, names, fontsize=8)
    axes[0].set_xlabel("Tissue macro-F1")
    axes[0].set_title("Tissue")
    _style_axes(axes[0])
    axes[1].barh(y, maes, color="#9a5b3b")
    axes[1].set_yticks(y, [""] * len(names))
    axes[1].set_xlabel("Age MAE (years)")
    axes[1].set_title("Age")
    _style_axes(axes[1])
    axes[2].barh(y, sexes, color="#3b9a6d")
    axes[2].set_yticks(y, [""] * len(names))
    axes[2].set_xlabel("Sex AUROC")
    axes[2].set_title("Sex")
    _style_axes(axes[2])
    fig.suptitle("Milestone 7G — methylation-input ranking (3 folds)", fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def write_analysis(
    path: Path,
    *,
    summary: dict[str, Any],
    means: list[dict[str, Any]],
    winner: dict[str, Any] | None,
) -> None:
    lines = [
        "# Stage 0 Milestone 7G — Methylation-only full evaluation",
        "",
        "Status: **done** when this report exists with Ranking + winner.",
        "",
        "## Budget and remaining ceilings",
        "",
        f"- Neural / classical loci: **{summary.get('max_loci')}** "
        f"(matrix has **{summary.get('n_loci_in_matrix')}**).",
        f"- Epochs: **{summary.get('max_epochs')}**; restarts: **{summary.get('n_restarts')}**.",
        "- Remaining ceiling: full **482 379** loci and a 2nd restart — deferred.",
        "- Trees: sklearn **HistGradientBoosting** (not LightGBM).",
        "- SVA: **10 train-fold PCA** surrogate variables (not Bioconductor `sva`).",
        "- Metadata-only omitted from ranking (7E′ leakage alarm only).",
        "",
        "## Topology under test",
        "",
        "```text",
        "CpG → typed region (gene | RBS) → RBS",
        "  ├─ allocated to gene → MBS",
        "  └─ orphan RBS",
        "leftover CpG (no typed region / former TBS) → direct",
        "late fusion: [orphan RBS | MBS | direct] → linear heads",
        "```",
        "",
        "No TBS arm (ADR 0009).",
        "",
    ]
    glossary_arms = [r["arm"] for r in means] + list(RANKING_ARMS)
    lines.extend(render_arm_glossary_section(glossary_arms))
    lines.extend(
        [
        "## Ranking (methylation-input only)",
        "",
        "Tissue macro-F1 / balanced accuracy exclude classes with zero training "
        "examples in that fold (study-grouped folds can hold an entire rare "
        "tissue class out of train; no model can predict those by "
        "construction, so counting them would penalize every arm for a fold-"
        "construction artifact rather than model quality). "
        "`n classes scored` / `n excluded (zero-shot)` make the denominator "
        "explicit per arm.",
        "",
        "| Arm | Tissue macro-F1 | Balanced acc | n classes scored | "
        "n excluded (zero-shot) | Sex AUROC | Age MAE | Age R² |",
        "|-----|-----------------|--------------|-------------------|"
        "------------------------|-----------|---------|--------|",
        ]
    )
    for r in means:
        lines.append(
            "| {arm} | {f1} | {bacc} | {nk} | {nex} | {sex} | {mae} | {r2} |".format(
                arm=r["arm"],
                f1=_fmt(r.get("tissue_macro_f1")),
                bacc=_fmt(r.get("tissue_balanced_accuracy")),
                nk=_fmt(r.get("tissue_n_classes_scored_mean")),
                nex=r.get("tissue_n_excluded_zero_shot_samples", 0),
                sex=_fmt(r.get("sex_auroc")),
                mae=_fmt(r.get("age_mae_years")),
                r2=_fmt(r.get("age_r2")),
            )
        )
    lines.extend(["", "## Winner (Milestone 7 topology)", ""])
    if winner:
        lines.append(
            f"**`{winner['arm']}`** — {winner.get('reason')}. "
            f"Tissue macro-F1={_fmt(winner.get('tissue_macro_f1'))}; "
            f"sex AUROC={_fmt(winner.get('sex_auroc'))}; "
            f"age MAE={_fmt(winner.get('age_mae_years'))}."
        )
    else:
        lines.append("_No ranking cells yet (cascade / classical incomplete)._")
    lines.extend(
        [
            "",
            "## ROC",
            "",
            "Sex AUROC and tissue one-vs-rest curves in `figures/` come from "
            "**neural fusion** scores (`N-cascade-l1`), not only HGB.",
            "",
            "## Artifacts",
            "",
            "- `summary.json` — full dump including sex metrics",
            "- `classical_baselines.json`",
            "- `transparent_baselines.json`",
            "- `arm_means.json`",
            "- `figures/roc_tissue_ovr_fusion.png`, `figures/roc_sex_fusion.png`",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def _fmt(x: float | None) -> str:
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "—"
    return f"{x:.3f}"


def load_cascade_summary(run_root: Path, report_summary: Path) -> dict[str, Any]:
    if report_summary.is_file():
        return json.loads(report_summary.read_text(encoding="utf-8"))
    folds = []
    for fold_dir in sorted(run_root.glob("fold_*")):
        metrics_path = fold_dir / "metrics.json"
        if metrics_path.is_file():
            folds.append(json.loads(metrics_path.read_text(encoding="utf-8")))
    return {"folds": folds, "arm": "N-cascade-l1"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "configs/experiment/stage0_7g_methylation_eval.yaml")
    parser.add_argument("--run-id", default=DEFAULT_RUN)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--skip-classical", action="store_true")
    parser.add_argument("--skip-transparent", action="store_true")
    parser.add_argument("--max-loci", type=int, default=None)
    args = parser.parse_args()

    paths = DataPaths.from_environment()
    cfg = load_experiment_config(
        args.config if args.config.is_absolute() else paths.project_root / args.config
    )
    max_loci = int(args.max_loci or cfg.get("cv_budget", {}).get("max_loci", MAX_LOCI))
    max_epochs = int(cfg.get("cv_budget", {}).get("max_epochs", 15))
    split_id = str(cfg.get("split_id", DEFAULT_SPLIT))
    report_dir = args.report_dir
    if not report_dir.is_absolute():
        report_dir = paths.project_root / report_dir
    report_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = report_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    folds_path = paths.artifact_root / "splits" / split_id / "folds.json"
    fold_pack = load_frozen_folds(folds_path)
    pheno_rel = Path(
        str(
            cfg.get("sample_phenotype_table")
            or "canonical/phenotypes/sample_phenotype_table_age_tissue_sex_full_v1.parquet"
        )
    )
    pheno_path = pheno_rel if pheno_rel.is_absolute() else paths.data_root / pheno_rel
    ont_path = paths.data_root / DEFAULT_ONTOLOGY
    _samples, phenotypes = samples_from_phenotype_table(pheno_path, ontology_path=ont_path)

    run_root = paths.artifact_root / "runs" / args.run_id
    cascade = load_cascade_summary(run_root, report_dir / "summary.json")
    if "max_loci" not in cascade:
        cascade["max_loci"] = max_loci
        cascade["max_epochs"] = max_epochs
        cascade["n_restarts"] = 1
        cascade["split_id"] = split_id

    classical_path = report_dir / "classical_baselines.json"
    if args.skip_classical and classical_path.is_file():
        classical = json.loads(classical_path.read_text(encoding="utf-8"))
    else:
        classical = run_classical_mvalue(
            data_root=paths.data_root,
            fold_pack=fold_pack,
            phenotypes=phenotypes,
            max_loci=max_loci,
        )
        write_json(classical_path, classical)

    transparent_path = report_dir / "transparent_baselines.json"
    if args.skip_transparent and transparent_path.is_file():
        transparent = json.loads(transparent_path.read_text(encoding="utf-8"))
    else:
        transparent = run_all_transparent_arms(
            data_root=paths.data_root,
            fold_pack=fold_pack,
            phenotypes=phenotypes,
            max_loci=max_loci,
        )
        write_json(transparent_path, transparent)

    rows = cascade_rows(cascade) + classical_rows(classical) + transparent_rows(transparent)
    means = aggregate_arms(rows)
    winner = pick_winner(means)
    write_json(report_dir / "arm_means.json", {"arms": means, "winner": winner})

    # Neural fusion ROC from fold 0 cascade metrics.
    fusion_roc_tissue = None
    fusion_sex = None
    if cascade.get("folds"):
        f0 = _metrics_from_fold_blob(cascade["folds"][0])
        fusion_roc_tissue = f0.get("tissue_roc")
        fusion_sex = f0.get("sex")
    plot_roc(
        fusion_roc_tissue,
        fusion_sex,
        fig_dir / "roc_tissue_ovr_fusion.png",
        fig_dir / "roc_sex_fusion.png",
        title_prefix="Neural fusion (N-cascade-l1), fold 0",
    )
    plot_ranking_bars(means, fig_dir / "ranking_bars.png")

    n_loci_matrix = classical.get("n_loci_in_matrix") or cascade.get("remaining_ceiling", {}).get(
        "n_loci_in_matrix"
    )
    summary = {
        "milestone": "7G",
        "topology": "rbs_gene_direct",
        "split_id": split_id,
        "run_id": args.run_id,
        "max_loci": max_loci,
        "max_epochs": max_epochs,
        "n_restarts": 1,
        "n_loci_in_matrix": n_loci_matrix,
        "remaining_ceiling": {
            "n_loci_used": max_loci,
            "n_loci_in_matrix": n_loci_matrix,
            "n_restarts_used": 1,
            "note": "Full-matrix (≈482k) train and 2nd restart deferred.",
        },
        "ranking": means,
        "winner": winner,
        "cascade": cascade,
        "metadata_only_in_ranking": False,
        "arms_present": [r["arm"] for r in means],
        "sex_metrics_present": any(r.get("sex_auroc") is not None for r in means),
    }
    write_json(report_dir / "summary.json", summary)
    write_analysis(report_dir / "analysis.md", summary=summary, means=means, winner=winner)
    print(f"wrote {report_dir}", flush=True)
    if winner:
        print(f"winner: {winner['arm']}", flush=True)


if __name__ == "__main__":
    main()
