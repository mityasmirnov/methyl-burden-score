#!/usr/bin/env python3
"""Write per-family + combined Hub real-matrix benchmark reports."""

from __future__ import annotations

import json
from pathlib import Path

FAMILIES = (
    ("age", "stage0-hub-age-studyholdout-v1", "matrix-hub-age-studyholdout-v1"),
    ("tissue", "stage0-hub-tissue-studyholdout-v1", "matrix-hub-tissue-studyholdout-v1"),
    ("blood", "stage0-hub-blood-studyholdout-v1", "matrix-hub-blood-studyholdout-v1"),
    ("brain", "stage0-hub-brain-studyholdout-v1", "matrix-hub-brain-studyholdout-v1"),
)


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _family_report(
    *,
    project_root: Path,
    family: str,
    run_id: str,
    matrix_id: str,
) -> dict:
    artifact_root = Path(
        __import__("os").environ.get(
            "MBS_ARTIFACT_ROOT", str(project_root / "artifacts")
        )
    )
    data_root = Path(
        __import__("os").environ.get("MBS_DATA_ROOT", str(project_root / "data"))
    )
    run_dir = artifact_root / "runs" / run_id
    matrix_dir = data_root / "canonical" / "matrices" / matrix_id
    metrics = _load_json(run_dir / "metrics.json")
    split = _load_json(run_dir / "split.json")
    study_subset = _load_json(matrix_dir / "study_subset.json")
    final = metrics.get("final") or {}
    external = metrics.get("external_test")
    payload = {
        "phenotype_family": family,
        "model_public_name": metrics.get("model_public_name", "deepMAT"),
        "run_id": run_id,
        "matrix_id": matrix_id,
        "matrix_path": str(matrix_dir),
        "study_ids": study_subset.get("study_ids"),
        "n_samples_matrix": study_subset.get("n_samples"),
        "split": {
            "split_id": split.get("split_id"),
            "mode": split.get("mode"),
            "train_studies": split.get("train_studies"),
            "validation_studies": split.get("validation_studies"),
            "external_test_studies": split.get("external_test_studies"),
            "n_train": len(split.get("train_sample_ids") or []),
            "n_validation": len(split.get("validation_sample_ids") or []),
            "n_external_test": len(split.get("external_test_sample_ids") or []),
        },
        "task": metrics.get("task"),
        "best_epoch": metrics.get("best_epoch"),
        "best_val_loss": metrics.get("best_val_loss"),
        "final_val_loss": final.get("val_loss"),
        "final_val_accuracy": final.get("val_accuracy"),
        "final_val_mae": final.get("val_mae"),
        "external_test": external,
        "tensorboard_dir": str(run_dir / "tb"),
        "metrics_jsonl": str(run_dir / "metrics.jsonl"),
        "synthetic_fixture": False,
    }
    return payload


def main() -> None:
    project_root = Path(__file__).resolve().parents[1]
    out_root = project_root / "reports" / "inspection" / "stage0_hub_real_benchmark"
    out_root.mkdir(parents=True, exist_ok=True)
    combined: list[dict] = []
    lines = [
        "# Stage 0 Hub real-matrix benchmark (study-grouped)",
        "",
        "Public model name: **deepMAT**. Package/CLI unchanged (`mbs`).",
        "",
        "These runs use **EWAS Data Hub profile-pack → canonical matrix** subsets,",
        "not synthetic fixtures. TensorBoard + `metrics.jsonl` are enabled on every run.",
        "",
    ]
    for family, run_id, matrix_id in FAMILIES:
        payload = _family_report(
            project_root=project_root,
            family=family,
            run_id=run_id,
            matrix_id=matrix_id,
        )
        combined.append(payload)
        fam_dir = out_root / family
        fam_dir.mkdir(parents=True, exist_ok=True)
        (fam_dir / "summary.json").write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        md = [
            f"# Hub real benchmark — {family}",
            "",
            f"- run_id: `{run_id}`",
            f"- matrix: `{matrix_id}`",
            f"- studies: `{', '.join(payload['study_ids'] or [])}`",
            f"- task: `{payload['task']}`",
            f"- split: train={payload['split']['train_studies']} "
            f"val={payload['split']['validation_studies']} "
            f"test={payload['split']['external_test_studies']}",
            f"- best_epoch: `{payload['best_epoch']}`",
            f"- best_val_loss: `{payload['best_val_loss']}`",
            f"- final_val_accuracy: `{payload['final_val_accuracy']}`",
            f"- final_val_mae: `{payload['final_val_mae']}`",
            f"- external_test: `{json.dumps(payload['external_test'])}`",
            "",
        ]
        (fam_dir / "summary.md").write_text("\n".join(md), encoding="utf-8")
        lines.extend(md)
        lines.append("")

    (out_root / "combined.json").write_text(
        json.dumps({"families": combined}, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    lines.extend(
        [
            "## Label harmonization / design notes",
            "",
            "- **blood** `phenotype_value` is per-sample cell-fraction strings — "
            "benchmark uses the `tissue` column instead.",
            "- **tissue / blood study-holdout:** each holdout study is a "
            "**single tissue type absent from train**. Zero holdout accuracy is "
            "expected for closed-set multiclass CE (unseen class IDs). This run "
            "validates leakage-safe splits + logging, not biological tissue "
            "prediction. Prefer multi-tissue studies or open-set / binary heads "
            "before interpreting biology.",
            "- **age** `final_val_mae` in run JSON is in **train-fold standardized** "
            "units; `external_test.mae` is destandardized years.",
            "- **disease** profile zip still downloading; age-pack GSM overlap for "
            "candidate studies is control-only (no case labels). Convert after pack "
            "completes; map empty→`control`; fix `ulcerative colitis` / "
            "`Ulcerative colitis` casing.",
            "- **cancer** profile zip incomplete — matrix convert deferred.",
            "",
        ]
    )
    (out_root / "summary.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {out_root}")


if __name__ == "__main__":
    main()
