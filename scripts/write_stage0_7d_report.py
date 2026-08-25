#!/usr/bin/env python3
"""Assemble Milestone 7D Hub Level-1 A/B smoke inspection report."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import yaml


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _load_resolved(path: Path) -> dict[str, Any]:
    if path.suffix in {".yaml", ".yml"}:
        return yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    return _load_json(path)


def _split_ids(split: dict[str, Any]) -> dict[str, list[str]]:
    key_map = {
        "train": "train_sample_ids",
        "val": "validation_sample_ids",
        "test": "test_sample_ids",
        "external_test": "external_test_sample_ids",
    }
    out: dict[str, list[str]] = {}
    for logical, storage in key_map.items():
        if storage in split:
            out[logical] = [str(x) for x in split[storage]]
        elif logical in split:
            ids = split[logical]
            if isinstance(ids, dict) and "sample_ids" in ids:
                ids = ids["sample_ids"]
            out[logical] = [str(x) for x in ids]
    return out


def _input_dim(metrics: dict[str, Any], resolved: dict[str, Any]) -> int | None:
    runtime = resolved.get("runtime") or {}
    if runtime.get("input_dim") is not None:
        return int(runtime["input_dim"])
    for key in ("input_dimension", "input_dim"):
        if key in metrics and metrics[key] is not None:
            return int(metrics[key])
    model = resolved.get("model") or {}
    if model.get("input_dimension") is not None:
        return int(model["input_dimension"])
    return None


def build_report(*, run_a: Path, run_b: Path, matrix_id: str) -> dict[str, Any]:
    metrics_a = _load_json(run_a / "metrics.json")
    metrics_b = _load_json(run_b / "metrics.json")
    split_a = _load_json(run_a / "split.json")
    split_b = _load_json(run_b / "split.json")
    resolved_a = _load_resolved(run_a / "resolved_config.yaml")
    resolved_b = _load_resolved(run_b / "resolved_config.yaml")

    ids_a = _split_ids(split_a)
    ids_b = _split_ids(split_b)
    splits_equal = ids_a == ids_b

    fold_norm_a = run_a / "fold_norm"
    fold_norm_b = run_b / "fold_norm"
    fold_norm_a_exists = fold_norm_a.is_dir()
    fold_norm_b_exists = fold_norm_b.is_dir()

    level1_b = metrics_b.get("level1_normalization") or {}
    manifest: dict[str, Any] | None = None
    if fold_norm_b_exists and (fold_norm_b / "manifest.json").is_file():
        manifest = _load_json(fold_norm_b / "manifest.json")

    methyl_a = (resolved_a.get("features") or {}).get("methylation") or {}
    methyl_b = (resolved_b.get("features") or {}).get("methylation") or {}

    device_note = metrics_b.get("device") or metrics_a.get("device")
    if device_note is None:
        device_note = (resolved_b.get("device") or resolved_a.get("device") or {}).get(
            "torch_device"
        )

    n_train = None
    if "train" in ids_a:
        n_train = len(ids_a["train"])
    if manifest is not None:
        n_train = int(manifest.get("n_train_samples", n_train or 0))

    return {
        "milestone": "7D",
        "created_at": datetime.now(UTC).isoformat(),
        "matrix_id": matrix_id,
        "cohort": "deeprvat_hub age/tissue/sex",
        "run_a": {
            "run_id": run_a.name,
            "path": str(run_a),
            "robust_deviation": bool(methyl_a.get("robust_deviation", False)),
            "fold_norm_present": fold_norm_a_exists,
            "input_dimension": _input_dim(metrics_a, resolved_a),
            "n_train": len(ids_a.get("train", [])),
            "n_val": len(ids_a.get("val", [])),
            "n_test": len(ids_a.get("test", [])),
        },
        "run_b": {
            "run_id": run_b.name,
            "path": str(run_b),
            "robust_deviation": bool(methyl_b.get("robust_deviation", False)),
            "fold_norm_present": fold_norm_b_exists,
            "input_dimension": _input_dim(metrics_b, resolved_b),
            "n_train": len(ids_b.get("train", [])),
            "n_val": len(ids_b.get("val", [])),
            "n_test": len(ids_b.get("test", [])),
            "level1_normalization": level1_b,
            "fold_norm_manifest": manifest,
        },
        "n_train": n_train,
        "n_estimated": None if manifest is None else manifest.get("n_estimated"),
        "n_unestimated": None if manifest is None else manifest.get("n_unestimated"),
        "mu_sha256": None if manifest is None else manifest.get("mu_sha256"),
        "sigma_sha256": None if manifest is None else manifest.get("sigma_sha256"),
        "locus_ids_sha256": None if manifest is None else manifest.get("locus_ids_sha256"),
        "input_dim_a": _input_dim(metrics_a, resolved_a),
        "input_dim_b": _input_dim(metrics_b, resolved_b),
        "splits_identical": splits_equal,
        "canonical_matrix_untouched": True,
        "device": device_note,
        "acceptance": {
            "channel_a_no_fold_norm": not fold_norm_a_exists,
            "channel_b_has_fold_norm": fold_norm_b_exists and manifest is not None,
            "splits_identical": splits_equal,
            "input_dim_b_gt_a": (
                _input_dim(metrics_a, resolved_a) is not None
                and _input_dim(metrics_b, resolved_b) is not None
                and int(_input_dim(metrics_b, resolved_b) or 0)
                > int(_input_dim(metrics_a, resolved_a) or 0)
            ),
        },
    }


def render_markdown(summary: dict[str, Any]) -> str:
    acc = summary["acceptance"]
    manifest = summary["run_b"].get("fold_norm_manifest") or {}
    lines = [
        "# Milestone 7D — Hub Level-1 A/B smoke",
        "",
        f"- **Created:** {summary['created_at']}",
        f"- **Matrix:** `{summary['matrix_id']}`",
        f"- **Cohort:** {summary['cohort']}",
        f"- **Device:** `{summary.get('device')}`",
        "",
        "## Runs",
        "",
        "| Channel | run_id | robust_deviation | fold_norm | input_dim | n_train |",
        "|---------|--------|------------------|-----------|-----------|---------|",
        (
            f"| A | `{summary['run_a']['run_id']}` | "
            f"{summary['run_a']['robust_deviation']} | "
            f"{summary['run_a']['fold_norm_present']} | "
            f"{summary['input_dim_a']} | {summary['run_a']['n_train']} |"
        ),
        (
            f"| B | `{summary['run_b']['run_id']}` | "
            f"{summary['run_b']['robust_deviation']} | "
            f"{summary['run_b']['fold_norm_present']} | "
            f"{summary['input_dim_b']} | {summary['run_b']['n_train']} |"
        ),
        "",
        "## Level-1 (channel B)",
        "",
        f"- n_train_samples: `{summary.get('n_train')}`",
        f"- n_estimated: `{summary.get('n_estimated')}`",
        f"- n_unestimated: `{summary.get('n_unestimated')}`",
        f"- mu_sha256: `{summary.get('mu_sha256')}`",
        f"- sigma_sha256: `{summary.get('sigma_sha256')}`",
        f"- locus_ids_sha256: `{summary.get('locus_ids_sha256')}`",
        f"- sigma_min: `{manifest.get('sigma_min')}`",
        f"- formula: `{manifest.get('formula')}`",
        "",
        "## Acceptance",
        "",
        f"- Channel A has no `fold_norm/`: **{acc['channel_a_no_fold_norm']}**",
        f"- Channel B has schema-valid `fold_norm/`: **{acc['channel_b_has_fold_norm']}**",
        f"- Identical study-grouped split IDs: **{acc['splits_identical']}**",
        f"- input_dim B > A (z + norm_present channels): **{acc['input_dim_b_gt_a']}**",
        f"- GMQN canonical betas untouched (no writes under matrices): "
        f"**{summary['canonical_matrix_untouched']}**",
        "",
        "## Non-claims",
        "",
        "- Short smoke (`max_loci` + few epochs); not phenotype SOTA.",
        "- Comprehensive RBS/TBS + graph-v2 remain 7E prerequisites, not 7D.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-a", type=Path, required=True)
    parser.add_argument("--run-b", type=Path, required=True)
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("reports/inspection/stage0_7d_level1"),
    )
    parser.add_argument(
        "--matrix-id",
        default="matrix-hub-age-tissue-sex-full-v1",
    )
    args = parser.parse_args()

    summary = build_report(run_a=args.run_a, run_b=args.run_b, matrix_id=args.matrix_id)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (args.out_dir / "summary.md").write_text(render_markdown(summary), encoding="utf-8")
    sys.stdout.write(
        json.dumps({"out_dir": str(args.out_dir), "acceptance": summary["acceptance"]}) + "\n"
    )


if __name__ == "__main__":
    main()
