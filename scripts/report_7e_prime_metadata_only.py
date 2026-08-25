#!/usr/bin/env python3
"""Write Milestone 7E′ metadata-only confounding ceilings (sidecar report)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from mbs.annotation.manifest import utc_now_iso, write_json
from mbs.paths import DataPaths
from mbs.training.controls import evaluate_metadata_only_ceiling
from mbs.training.loop import maybe_constrained_split
from mbs.training.phenotype_table import load_tissue_ontology
from mbs.training.phenotypes import load_longform_multilabel, load_multitask_phenotypes


def _pheno_by_ids(
    phenotypes: list[Any],
    ids: set[str],
) -> list[Any]:
    return [p for p in phenotypes if p.sample_id in ids]


def _ceiling_for_table(
    *,
    label: str,
    phenotype_table: Path,
    tissue_ontology: Path | None,
    split: dict[str, Any],
    disease_matrix: Path | None = None,
    cancer_matrix: Path | None = None,
) -> dict[str, Any]:
    class_names = None
    if tissue_ontology is not None and tissue_ontology.is_file():
        class_names = load_tissue_ontology(tissue_ontology).class_names
    phenotypes, _ = load_multitask_phenotypes(
        phenotype_table,
        class_names=class_names,
    )
    train_ids = set(split["train_sample_ids"])
    val_ids = set(split.get("validation_sample_ids") or [])
    test_ids = set(split.get("external_test_sample_ids") or [])
    train = _pheno_by_ids(phenotypes, train_ids)
    eval_sets: dict[str, list[Any]] = {}
    if val_ids:
        eval_sets["validation"] = _pheno_by_ids(phenotypes, val_ids)
    if test_ids:
        eval_sets["external_test"] = _pheno_by_ids(phenotypes, test_ids)

    disease_maps = None
    cancer_maps = None
    all_ids = sorted(train_ids | val_ids | test_ids)
    if disease_matrix is not None and disease_matrix.is_file():
        disease_maps = load_longform_multilabel(
            disease_matrix,
            sample_ids=all_ids,
            min_count=5,
        )
    if cancer_matrix is not None and cancer_matrix.is_file():
        cancer_maps = load_longform_multilabel(
            cancer_matrix,
            sample_ids=all_ids,
            min_count=5,
        )

    metrics = evaluate_metadata_only_ceiling(
        train=train,
        eval_sets=eval_sets,
        disease_maps=disease_maps,
        cancer_maps=cancer_maps,
    )
    return {
        "label": label,
        "split_id": split.get("split_id"),
        "n_train": len(train),
        "n_validation": len(eval_sets.get("validation", [])),
        "n_external_test": len(eval_sets.get("external_test", [])),
        "metrics": metrics,
    }


def main() -> None:
    paths = DataPaths.from_environment()
    report_dir = paths.project_root / "reports" / "inspection" / "stage0_7e_prime"
    report_dir.mkdir(parents=True, exist_ok=True)

    ats_table = (
        paths.data_root
        / "canonical"
        / "phenotypes"
        / "sample_phenotype_table_age_tissue_sex_full_v1.parquet"
    )
    ats_ont = (
        paths.data_root
        / "canonical"
        / "phenotypes"
        / "tissue_ontology_age_tissue_sex_full_v1.yaml"
    )
    ats_split_path = (
        paths.artifact_root / "runs" / "stage0-flat-deeprvat-age-tissue-sex-full-v1" / "split.json"
    )
    ats_split = json.loads(ats_split_path.read_text(encoding="utf-8"))

    hub_table = (
        paths.data_root
        / "canonical"
        / "phenotypes"
        / "sample_phenotype_table_hub_nine_pack_v1.parquet"
    )
    hub_ont = (
        paths.data_root / "canonical" / "phenotypes" / "tissue_ontology_hub_nine_pack_v1.yaml"
    )
    hub_frame = pd.read_parquet(hub_table)
    sample_rows = [
        {
            "sample_id": str(r["sample_id"]),
            "study_id": str(r["study_id"]),
            "platform": r.get("platform"),
        }
        for r in hub_frame.to_dict(orient="records")
    ]
    hub_split = maybe_constrained_split(
        sample_rows,
        seed=42,
        train_fraction=0.7,
        val_fraction=0.15,
        split_id="hub-nine-pack-full-auto-v1",
    )
    (report_dir / "hub_nine_pack_split.json").write_text(
        json.dumps(hub_split, indent=2) + "\n",
        encoding="utf-8",
    )

    matrices = paths.data_root / "canonical" / "matrices"
    disease_pheno = matrices / "matrix-hub-disease-full-v1" / "sample_phenotypes.parquet"
    cancer_pheno = matrices / "matrix-hub-cancer-full-v1" / "sample_phenotypes.parquet"

    ats_report = _ceiling_for_table(
        label="ats_freeze_hub-age-tissue-sex-full-auto-v1",
        phenotype_table=ats_table,
        tissue_ontology=ats_ont,
        split=ats_split,
    )
    hub_report = _ceiling_for_table(
        label="hub_nine_pack_full-auto-v1",
        phenotype_table=hub_table,
        tissue_ontology=hub_ont,
        split=hub_split,
        disease_matrix=disease_pheno,
        cancer_matrix=cancer_pheno,
    )

    payload = {
        "generated_at": utc_now_iso(),
        "protocol": "fit_train_score_holdout",
        "note": (
            "Metadata-only confounding ceiling (study/platform/tissue one-hots). "
            "Sidecar for 7E-prime; 7E bake-off report must include the ATS ceiling. "
            "Not a replacement for neural training."
        ),
        "cohorts": [ats_report, hub_report],
        "artifacts": {
            "virtual_matrix_id": "matrix-hub-nine-pack-virtual-v1",
            "hub_phenotype_table": str(hub_table),
            "hub_split_id": "hub-nine-pack-full-auto-v1",
        },
    }
    write_json(report_dir / "metadata_only.json", payload)

    lines = [
        "# Milestone 7E-prime metadata-only ceiling",
        "",
        f"- Generated: `{payload['generated_at']}`",
        f"- Protocol: `{payload['protocol']}`",
        "",
        payload["note"],
        "",
    ]
    for cohort in payload["cohorts"]:
        lines.extend(
            [
                f"## {cohort['label']}",
                "",
                f"- split_id: `{cohort['split_id']}`",
                f"- n_train / val / test: **{cohort['n_train']}** / "
                f"**{cohort['n_validation']}** / **{cohort['n_external_test']}**",
                "",
            ]
        )
        metrics = cohort["metrics"]
        for fold_name in ("validation", "external_test"):
            fold = metrics.get(fold_name)
            if not isinstance(fold, dict):
                continue
            lines.append(f"### {fold_name}")
            lines.append("")
            for task, vals in fold.items():
                if not isinstance(vals, dict):
                    continue
                summary = ", ".join(f"{k}={v}" for k, v in vals.items() if not isinstance(v, dict))
                lines.append(f"- **{task}**: {summary}")
            lines.append("")
    lines.extend(
        [
            "## Virtual Hub store",
            "",
            "- `matrix-hub-nine-pack-virtual-v1` (route + indices; no dense Zarr)",
            "- Phenotype table: `sample_phenotype_table_hub_nine_pack_v1.parquet`",
            "- Blood `cell_component` is **not** a pack-wide head",
            "",
        ]
    )
    (report_dir / "metadata_only.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    # Operator-facing path for shell pipelines.
    print(json.dumps({"report_dir": str(report_dir)}, indent=2))  # noqa: T201


if __name__ == "__main__":
    main()
