#!/usr/bin/env python3
"""Write Milestone 6 hierarchical vs flat 5d inspection report."""

from __future__ import annotations

import json
import os
from pathlib import Path

DEFAULT_RUN_ID = "stage0-hier-deeprvat-age-tissue-sex-full-v1"
FALLBACK_RUN_ID = "stage0-hier-smoke-maxloci"
FLAT_RUN_ID = "stage0-flat-deeprvat-age-tissue-sex-full-v1"
MATRIX_ID = "matrix-hub-age-tissue-sex-full-v1"


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _safe_metrics(path: Path) -> dict | None:
    if not path.is_file():
        return None
    return _load_json(path)


def main() -> None:
    project_root = Path(__file__).resolve().parents[1]
    artifact_root = Path(os.environ.get("MBS_ARTIFACT_ROOT", str(project_root / "artifacts")))
    preferred = os.environ.get("MBS_HIER_RUN_ID", DEFAULT_RUN_ID)
    run_id = preferred
    run_dir = artifact_root / "runs" / run_id
    metrics = _safe_metrics(run_dir / "metrics.json")
    if metrics is None or metrics.get("overfit_fixture"):
        run_id = FALLBACK_RUN_ID
        run_dir = artifact_root / "runs" / run_id
        metrics = _safe_metrics(run_dir / "metrics.json")
    ckpt_dir = artifact_root / "checkpoints" / run_id
    flat_run = artifact_root / "runs" / FLAT_RUN_ID
    split = _safe_metrics(run_dir / "split.json")
    flat_metrics = _safe_metrics(flat_run / "metrics.json")
    flat_report = _safe_metrics(
        project_root / "reports" / "inspection" / "stage0_5d_max_n" / "summary.json"
    )

    if metrics is None:
        raise SystemExit(
            f"missing hierarchical metrics under {artifact_root / 'runs'}; "
            "run: uv run mbs train hierarchical "
            "--config configs/experiment/stage0_hier_deeprvat_full.yaml"
        )

    final = metrics.get("final") or {}
    external = metrics.get("external_test") or {}
    ablations = metrics.get("ablations") or {}
    flat_external = (flat_metrics or {}).get("external_test") or {}
    if not flat_external and flat_report:
        flat_external = flat_report.get("external_test") or {}

    payload = {
        "milestone": "6",
        "model_public_name": metrics.get("model_public_name", "deepMAT-hierarchical"),
        "run_id": run_id,
        "preferred_run_id": preferred,
        "config": "configs/experiment/stage0_hier_deeprvat_full.yaml",
        "plan": "docs/plans/milestone-6-hierarchical-region-model.md",
        "matrix_id": MATRIX_ID,
        "flat_baseline_run_id": FLAT_RUN_ID,
        "max_loci": (split or {}).get("max_loci"),
        "split": {
            "split_id": (split or {}).get("split_id"),
            "reused_flat_split": (split or {}).get("reused_flat_split"),
            "flat_split_run_id": (split or {}).get("flat_split_run_id"),
            "n_train": len((split or {}).get("train_sample_ids") or []),
            "n_validation": len((split or {}).get("validation_sample_ids") or []),
            "n_external_test": len((split or {}).get("external_test_sample_ids") or []),
        },
        "topology": {
            "n_genes": metrics.get("n_genes"),
            "n_regions": metrics.get("n_regions"),
            "n_typed_edges": metrics.get("n_typed_edges"),
            "n_unassigned_regions": metrics.get("n_unassigned_regions"),
            "region_types": metrics.get("region_types"),
            "unassigned_policy": (
                "singleton region_type=unassigned → synthetic gene __unassigned__; "
                "Illumina-coordinate-unmapped probes remain matrix-excluded"
            ),
        },
        "task": metrics.get("task"),
        "n_classes": metrics.get("n_classes"),
        "best_epoch": metrics.get("best_epoch"),
        "best_val_loss": metrics.get("best_val_loss"),
        "final": final,
        "external_test": external,
        "ablations": ablations,
        "comparison_to_flat_5d": {
            "flat_external_test": flat_external,
            "hier_external_test": external,
            "note": (
                "Same preferred split when reused_flat_split=true; "
                "compare tissue accuracy / age MAE / sex accuracy"
            ),
        },
        "checkpoints": {
            "best": str(ckpt_dir / "best.pt"),
            "last": str(ckpt_dir / "last.pt"),
            "manifest": str(ckpt_dir / "checkpoint_manifest.json"),
        },
        "tensorboard_dir": str(run_dir / "tb"),
        "metrics_jsonl": str(run_dir / "metrics.jsonl"),
    }

    out_dir = project_root / "reports" / "inspection" / "stage0_6_hierarchical"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "summary.json").write_text(
        json.dumps(payload, indent=2, default=str) + "\n",
        encoding="utf-8",
    )

    flat_acc = flat_external.get("accuracy")
    hier_acc = external.get("accuracy")
    flat_mae = flat_external.get("mae")
    hier_mae = external.get("mae")
    max_loci = payload.get("max_loci")
    lines = [
        "# Milestone 6 — hierarchical region model",
        "",
        f"- **run_id:** `{run_id}`",
        f"- **config:** `{payload['config']}`",
        f"- **matrix:** `{MATRIX_ID}`",
        f"- **max_loci:** `{max_loci}` (null = full matrix columns)",
        f"- **plan:** `{payload['plan']}`",
        "",
        "## Topology",
        "",
        f"- genes: `{payload['topology']['n_genes']}`",
        f"- regions: `{payload['topology']['n_regions']}`",
        f"- typed edges: `{payload['topology']['n_typed_edges']}`",
        f"- unassigned singleton regions: `{payload['topology']['n_unassigned_regions']}`",
        f"- region types: `{payload['topology']['region_types']}`",
        f"- policy: {payload['topology']['unassigned_policy']}",
        "",
        "## Split",
        "",
        f"- reused flat 5d split: `{payload['split']['reused_flat_split']}`",
        f"- train/val/test n: "
        f"{payload['split']['n_train']}/"
        f"{payload['split']['n_validation']}/"
        f"{payload['split']['n_external_test']}",
        "",
        "## External / holdout metrics",
        "",
        f"- hierarchical: accuracy={hier_acc}, mae={hier_mae}, "
        f"sex_accuracy={external.get('sex_accuracy')}",
        f"- flat 5d: accuracy={flat_acc}, mae={flat_mae}, "
        f"sex_accuracy={flat_external.get('sex_accuracy')}",
        "",
        "## Ablations (holdout subset)",
        "",
    ]
    for name, row in ablations.items():
        lines.append(
            f"- **{name}:** accuracy={row.get('accuracy')}, mae={row.get('mae')}, "
            f"sex_accuracy={row.get('sex_accuracy')}, n={row.get('n_samples')}"
        )
    lines.extend(
        [
            "",
            "## Artifacts",
            "",
            f"- metrics: `{run_dir / 'metrics.json'}`",
            f"- checkpoints: `{ckpt_dir}`",
            f"- TensorBoard: `{run_dir / 'tb'}`",
            "- full uncapped train log: `scratch/logs/hier_full.log`",
            "",
        ]
    )
    (out_dir / "summary.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {out_dir / 'summary.md'} (run_id={run_id})")  # noqa: T201


if __name__ == "__main__":
    main()
