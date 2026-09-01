#!/usr/bin/env python3
"""Aggregate 7G cascade tissue probe arms into inspection report."""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from mbs.annotation.manifest import write_json
from mbs.paths import DataPaths
from mbs.training.dev_cv import load_frozen_folds, samples_from_phenotype_table
from mbs.training.loop import load_experiment_config

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "configs/experiment/stage0_7g_cascade_tissue_probe.yaml"
DEFAULT_REPORT = ROOT / "reports/inspection/stage0_7g_cascade_tissue_probe"
DEFAULT_ONTOLOGY = "canonical/phenotypes/tissue_ontology_age_tissue_sex_full_v1.yaml"
F1_GATE = 0.20

ARM_LABELS = {
    "P0-baseline": "Replay 7G N-cascade-l1 (late fusion)",
    "P1-fusion-tissue-heavy": "Balanced logistic + PCA(32) on saved scores",
    "P2-end2end-tissue-weight": "tissue_loss_weight=3, age_loss_weight=0.3",
    "P3-region-head-bypass": "T-mean-region transparent (no cascade train)",
    "P2-fusion-balanced": "P2 scores + balanced logistic fusion",
    "P4-pooling-mean": "P2 weights + mean/mean pooling",
    "P4-fusion-balanced": "P4 scores + balanced logistic fusion",
    "P5-epochs-30": "P2 weights + 30-epoch cap + early stop",
    "P5-fusion-balanced": "P5 scores + balanced logistic fusion",
}


def _mean_std(xs: list[float | None]) -> tuple[float | None, float | None]:
    vals = [float(x) for x in xs if x is not None and not (isinstance(x, float) and math.isnan(x))]
    if not vals:
        return None, None
    arr = np.asarray(vals, dtype=np.float64)
    return float(arr.mean()), float(arr.std(ddof=1)) if arr.size > 1 else 0.0


def _fmt(x: float | None) -> str:
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "—"
    return f"{x:.3f}"


def _metrics_from_fold_blob(blob: dict[str, Any]) -> dict[str, Any]:
    m = blob.get("metrics") if isinstance(blob.get("metrics"), dict) else blob
    age = m.get("age") if isinstance(m.get("age"), dict) else None
    tissue = m.get("tissue") if isinstance(m.get("tissue"), dict) else None
    sex = m.get("sex") if isinstance(m.get("sex"), dict) else None
    return {
        "tissue_macro_f1": (tissue or {}).get("macro_f1"),
        "tissue_balanced_accuracy": (tissue or {}).get("balanced_accuracy"),
        "age_mae": (age or {}).get("mae"),
        "age_r2": (age or {}).get("r2"),
        "sex_auroc": (sex or {}).get("auroc"),
        "checkpoint_selection": blob.get("checkpoint_selection"),
    }


def aggregate_arm(arm_id: str, folds: list[dict[str, Any]]) -> dict[str, Any]:
    rows = [_metrics_from_fold_blob(f) for f in folds]
    f1_mean, f1_std = _mean_std([r["tissue_macro_f1"] for r in rows])
    bacc_mean, _ = _mean_std([r["tissue_balanced_accuracy"] for r in rows])
    mae_mean, _ = _mean_std([r["age_mae"] for r in rows])
    r2_mean, _ = _mean_std([r["age_r2"] for r in rows])
    sex_mean, _ = _mean_std([r["sex_auroc"] for r in rows])
    per_fold = [
        {
            "fold_id": folds[i].get("fold_id", i) if i < len(folds) else i,
            **rows[i],
            "checkpoint_selection": folds[i].get("checkpoint_selection"),
        }
        for i in range(len(rows))
    ]
    return {
        "arm_id": arm_id,
        "label": ARM_LABELS.get(arm_id, arm_id),
        "tissue_macro_f1": f1_mean,
        "tissue_macro_f1_std": f1_std,
        "tissue_balanced_accuracy": bacc_mean,
        "age_mae_years": mae_mean,
        "age_r2": r2_mean,
        "sex_auroc": sex_mean,
        "per_fold": per_fold,
        "n_folds": len(folds),
    }


def load_locked_tissue_comparator(report_dir: Path) -> dict[str, Any] | None:
    """Load the committed 7G C-mvalue-enet cells for paired context."""
    path = report_dir.parent / "stage0_7g_methylation_eval" / "classical_baselines.json"
    if not path.is_file():
        return None
    payload = json.loads(path.read_text(encoding="utf-8"))
    per_fold: list[dict[str, Any]] = []
    for fold in payload.get("folds") or []:
        blob = (fold.get("arms") or {}).get("C-mvalue-enet") or {}
        tissue = blob.get("tissue") or {}
        per_fold.append(
            {
                "fold_id": fold.get("fold"),
                "tissue_macro_f1": tissue.get("macro_f1"),
            }
        )
    mean, std = _mean_std([row.get("tissue_macro_f1") for row in per_fold])
    return {
        "arm": "C-mvalue-enet",
        "tissue_macro_f1": mean,
        "tissue_macro_f1_std": std,
        "per_fold": per_fold,
    }


def study_composition_table(
    fold_pack: dict[str, Any],
    phenotypes: list[Any],
) -> list[dict[str, Any]]:
    """Per-fold train/test tissue class counts (read-only confounding check)."""
    ph_by_id = {p.sample_id: p for p in phenotypes}
    rows: list[dict[str, Any]] = []
    for fold_i, fold in enumerate(fold_pack["folds"]):
        for split_name, ids in (
            ("train", fold.get("train_sample_ids") or []),
            ("test", fold.get("external_test_sample_ids") or fold.get("validation_sample_ids") or []),
        ):
            tissues = [
                int(ph_by_id[s].class_index)
                for s in ids
                if s in ph_by_id and ph_by_id[s].tissue_mask
            ]
            studies = [
                str(ph_by_id[s].study_id or "NA")
                for s in ids
                if s in ph_by_id and ph_by_id[s].tissue_mask
            ]
            rows.append(
                {
                    "fold_id": fold.get("fold_id", fold_i),
                    "split": split_name,
                    "n_tissue_labeled": len(tissues),
                    "n_unique_studies": len(set(studies)),
                    "top_tissue_counts": dict(Counter(tissues).most_common(8)),
                    "top_study_counts": dict(Counter(studies).most_common(8)),
                }
            )
    return rows


def diagnose(arm_means: dict[str, dict[str, Any]]) -> tuple[str, list[str]]:
    """Return primary failure mode and bullet evidence."""
    p0 = arm_means.get("P0-baseline", {}).get("tissue_macro_f1")
    p1 = arm_means.get("P1-fusion-tissue-heavy", {}).get("tissue_macro_f1")
    p2 = arm_means.get("P2-end2end-tissue-weight", {}).get("tissue_macro_f1")
    p3 = arm_means.get("P3-region-head-bypass", {}).get("tissue_macro_f1")
    bullets: list[str] = []
    if p0 is not None:
        bullets.append(f"P0 baseline tissue macro-F1 = {_fmt(p0)} (7G replay).")
    if p1 is not None and p0 is not None and p1 - p0 >= 0.05:
        bullets.append(
            f"P1 fusion-heavy ({_fmt(p1)}) lifts tissue F1 >=0.05 vs P0 -> **fusion bottleneck**."
        )
        return "fusion_bottleneck", bullets
    if p2 is not None and p0 is not None and p2 - p0 >= 0.05:
        bullets.append(
            f"P2 reweighted training ({_fmt(p2)}) lifts tissue F1 >=0.05 vs P0 "
            "-> **task competition**."
        )
        return "task_competition", bullets
    if p3 is not None and p3 >= 0.25 and (p0 is None or p0 < 0.15):
        bullets.append(
            f"P3 region-mean ({_fmt(p3)}) near classical ~0.33 while cascade arms stay low -> "
            "**gene aggregation (max-pool MBS) loses tissue signal**."
        )
        return "aggregation_loss", bullets
    bullets.append(
        "No single arm cleared the planned lift thresholds; likely mixed capacity / confounding."
    )
    return "inconclusive", bullets


def milestone7_recommendation(
    arm_means: dict[str, dict[str, Any]],
    diagnosis: str,
) -> str:
    best_f1 = max(
        (v.get("tissue_macro_f1") or -1.0 for v in arm_means.values()),
        default=-1.0,
    )
    if best_f1 >= F1_GATE and diagnosis in {"fusion_bottleneck", "task_competition"}:
        return (
            "**Proceed with narrowed hyperparameter grid** (P4-P5 + fusion solver sweep) "
            "before Milestone 7 OOF; cascade may be salvageable for tissue with fusion/loss fixes."
        )
    if best_f1 < F1_GATE:
        return (
            f"**Draft ADR before full OOF spend**: mean probe tissue F1 ({_fmt(best_f1)}) "
            "stays below ~0.20. Export 7F product scores (MBS / orphan RBS / direct) in "
            "Milestone 7 while tissue **benchmarking** may remain `C-mvalue-enet`."
        )
    return (
        "**Proceed cascade OOF as-is** for product score export; keep classical enet as "
        "methylation comparator winner for tissue phenotype tables."
    )


def checkpoint_audit(arm_means: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for arm_id in (
        "P0-baseline",
        "P2-end2end-tissue-weight",
        "P4-pooling-mean",
        "P5-epochs-30",
    ):
        arm = arm_means.get(arm_id)
        if not arm:
            continue
        for fold in arm.get("per_fold") or []:
            sel = fold.get("checkpoint_selection") or {}
            rows.append(
                {
                    "arm_id": arm_id,
                    "fold_id": fold.get("fold_id"),
                    "best_epoch": sel.get("best_epoch"),
                    "selection": sel.get("selection"),
                    "val_tissue_macro_f1_at_best": (
                        (sel.get("best_validation_metrics") or {}).get("tissue_macro_f1")
                    ),
                    "epochs_completed": sel.get("epochs_completed"),
                    "early_stopping": sel.get("early_stopping"),
                }
            )
    return rows


def plot_tissue_bars(arm_means: list[dict[str, Any]], path: Path) -> None:
    names = [ARM_LABELS.get(a["arm_id"], a["arm_id"]) for a in arm_means]
    f1s = [a.get("tissue_macro_f1") or float("nan") for a in arm_means]
    stds = [a.get("tissue_macro_f1_std") or 0.0 for a in arm_means]
    y = np.arange(len(names))
    fig, ax = plt.subplots(figsize=(8.5, max(4.5, 0.55 * len(names))))
    ax.barh(y, f1s, xerr=stds, color="#3b6d9a", capsize=3)
    ax.set_yticks(y, names, fontsize=9)
    ax.set_xlabel("Tissue macro-F1 (mean ± std, 3 folds)")
    ax.axvline(F1_GATE, color="#c44", linestyle="--", linewidth=1, label=f"gate {F1_GATE}")
    ax.legend(loc="lower right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", linestyle=":", alpha=0.4)
    fig.suptitle("7G cascade tissue probe", fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def write_analysis(
    path: Path,
    *,
    arm_means: list[dict[str, Any]],
    diagnosis: str,
    diagnosis_bullets: list[str],
    recommendation: str,
    study_rows: list[dict[str, Any]],
    ckpt_rows: list[dict[str, Any]],
    locked_comparator: dict[str, Any] | None,
    cfg: dict[str, Any],
) -> None:
    lines = [
        "# Stage 0 - 7G cascade tissue probe (P0-P5)",
        "",
        "Frozen split `hub-ats-7e-3fold-v1`; **65536 loci / 1 restart**. "
        "P0-P4 use 15 epochs; P5 uses a 30-epoch ceiling with validation-tissue early stopping.",
        "",
        "## Per-arm summary",
        "",
        "| Arm | Tissue macro-F1 | Balanced acc | Sex AUROC | Age MAE | Age R² |",
        "|-----|-----------------|--------------|-----------|---------|--------|",
    ]
    for a in arm_means:
        lines.append(
            f"| {a['arm_id']} | {_fmt(a.get('tissue_macro_f1'))} | "
            f"{_fmt(a.get('tissue_balanced_accuracy'))} | {_fmt(a.get('sex_auroc'))} | "
            f"{_fmt(a.get('age_mae_years'))} | {_fmt(a.get('age_r2'))} |"
        )
    if locked_comparator is not None:
        lines.extend(
            [
                "",
                "## Locked 7G tissue comparator",
                "",
                f"`C-mvalue-enet`: mean F1 "
                f"{_fmt(locked_comparator.get('tissue_macro_f1'))}. "
                "It remains the locked phenotype comparator until the 7H same-panel benchmark; "
                "post-7G targeted arms are reported as salvage evidence, not retroactively "
                "inserted into the 7G winner selection.",
            ]
        )
    lines.extend(["", "## Per-fold tissue macro-F1", ""])
    for a in arm_means:
        lines.append(f"### {a['arm_id']}")
        for fold in a.get("per_fold") or []:
            lines.append(
                f"- fold {fold.get('fold_id')}: F1={_fmt(fold.get('tissue_macro_f1'))}, "
                f"sex AUROC={_fmt(fold.get('sex_auroc'))}, age MAE={_fmt(fold.get('age_mae'))}"
            )
        lines.append("")
    lines.extend(["## Diagnosis", "", f"**Primary:** `{diagnosis}`", ""])
    lines.extend(diagnosis_bullets)
    lines.extend(["", "## Checkpoint audit (trained cascade arms)", ""])
    if ckpt_rows:
        lines.append(
            "| Arm | Fold | Best epoch | Epochs run | Best val tissue F1 | "
            "Selection | Early stopped |"
        )
        lines.append(
            "|-----|------|------------|------------|---------------------|"
            "-----------|---------------|"
        )
        for r in ckpt_rows:
            lines.append(
                f"| {r['arm_id']} | {r['fold_id']} | {r.get('best_epoch', '—')} | "
                f"{r.get('epochs_completed', '—')} | "
                f"{_fmt(r.get('val_tissue_macro_f1_at_best'))} | "
                f"{r.get('selection', '—')} | "
                f"{(r.get('early_stopping') or {}).get('stopped_early', False)} |"
            )
    else:
        lines.append("_No checkpoint metadata found._")
    lines.extend(["", "## Study composition (tissue-labeled)", ""])
    for r in study_rows:
        lines.append(
            f"- fold {r['fold_id']} {r['split']}: {r['n_tissue_labeled']} samples, "
            f"{r['n_unique_studies']} studies; top tissues: {r['top_tissue_counts']}"
        )
    lines.extend(
        [
            "",
            "## Milestone 7 recommendation",
            "",
            recommendation,
            "",
            "## Product scores versus phenotype comparator",
            "",
            "The product export remains the 7F cascade: **MBS + qualified orphan RBS + "
            "direct loci**. Tissue ranking is a separate phenotype benchmark; "
            "`C-mvalue-enet` remains the locked classical comparator until a cascade "
            "arm beats it under the same folds and input panel. See ADR 0010.",
            "",
            "## Artifacts",
            "",
            f"- Config: `{cfg.get('experiment', {}).get('name', 'stage0_7g_cascade_tissue_probe')}`",
            "- `arm_means.json`, `per_arm/*.json`, `figures/tissue_f1_bars.png`",
            "",
            "Phase-2 is complete only when P4 and P5 each have all three fold artifacts. "
            "Do not infer their performance from P0-P3.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    paths = DataPaths.from_environment()
    config_path = args.config if args.config.is_absolute() else paths.project_root / args.config
    cfg = load_experiment_config(config_path)
    report_dir = args.report_dir
    if not report_dir.is_absolute():
        report_dir = paths.project_root / report_dir
    per_arm_dir = report_dir / "per_arm"
    fig_dir = report_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    split_id = str(cfg.get("split_id", "hub-ats-7e-3fold-v1"))
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

    arm_means_list: list[dict[str, Any]] = []
    arm_means_map: dict[str, dict[str, Any]] = {}
    for path in sorted(per_arm_dir.glob("*.json")):
        payload = json.loads(path.read_text(encoding="utf-8"))
        arm_id = str(payload.get("arm_id", path.stem))
        folds = payload.get("folds") or []
        agg = aggregate_arm(arm_id, folds)
        arm_means_list.append(agg)
        arm_means_map[arm_id] = agg

    diagnosis, bullets = diagnose(arm_means_map)
    recommendation = milestone7_recommendation(arm_means_map, diagnosis)
    study_rows = study_composition_table(fold_pack, phenotypes)
    ckpt_rows = checkpoint_audit(arm_means_map)
    locked_comparator = load_locked_tissue_comparator(report_dir)

    write_json(
        report_dir / "arm_means.json",
        {
            "arms": arm_means_list,
            "diagnosis": diagnosis,
            "recommendation": recommendation,
            "f1_gate": F1_GATE,
            "study_composition": study_rows,
            "checkpoint_audit": ckpt_rows,
            "locked_tissue_comparator": locked_comparator,
        },
    )
    plot_tissue_bars(arm_means_list, fig_dir / "tissue_f1_bars.png")
    write_analysis(
        report_dir / "analysis.md",
        arm_means=arm_means_list,
        diagnosis=diagnosis,
        diagnosis_bullets=bullets,
        recommendation=recommendation,
        study_rows=study_rows,
        ckpt_rows=ckpt_rows,
        locked_comparator=locked_comparator,
        cfg=cfg,
    )
    print(f"wrote {report_dir / 'analysis.md'}", flush=True)


if __name__ == "__main__":
    main()
