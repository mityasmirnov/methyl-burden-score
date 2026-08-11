#!/usr/bin/env python3
"""Write Milestone 5d max-N DeepRVAT flat baseline inspection report."""

from __future__ import annotations

import json
from pathlib import Path

import pyarrow.parquet as pq

RUN_ID = "stage0-flat-deeprvat-age-tissue-sex-full-v1"
MATRIX_ID = "matrix-hub-age-tissue-sex-full-v1"
SOURCE_MATRICES = (
    "matrix-hub-age-full-v1",
    "matrix-hub-tissue-full-v1",
    "matrix-hub-sex-full-v1",
)


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _matrix_shape(data_root: Path, matrix_id: str) -> dict:
    man = _load_json(data_root / "canonical" / "matrices" / matrix_id / "matrix_manifest.json")
    shape = man.get("shape") or [None, None]
    return {
        "matrix_id": matrix_id,
        "shape": shape,
        "n_samples": shape[0],
        "n_loci": shape[1],
        "platform_id": man.get("platform_id"),
        "genome_build": man.get("genome_build"),
        "created_at": man.get("created_at"),
        "notes": man.get("notes"),
    }


def main() -> None:
    project_root = Path(__file__).resolve().parents[1]
    artifact_root = Path(
        __import__("os").environ.get("MBS_ARTIFACT_ROOT", str(project_root / "artifacts"))
    )
    data_root = Path(__import__("os").environ.get("MBS_DATA_ROOT", str(project_root / "data")))
    run_dir = artifact_root / "runs" / RUN_ID
    ckpt_dir = artifact_root / "checkpoints" / RUN_ID
    metrics = _load_json(run_dir / "metrics.json")
    split = _load_json(run_dir / "split.json")
    ph_path = (
        data_root
        / "canonical"
        / "phenotypes"
        / "sample_phenotype_table_age_tissue_sex_full_v1.parquet"
    )
    table = pq.read_table(ph_path)
    df = table.to_pandas()
    sources = [_matrix_shape(data_root, mid) for mid in SOURCE_MATRICES]
    merged = _matrix_shape(data_root, MATRIX_ID)
    final = metrics.get("final") or {}
    external = metrics.get("external_test") or {}
    payload = {
        "milestone": "5d",
        "model_public_name": metrics.get("model_public_name", "deepMAT"),
        "run_id": RUN_ID,
        "config": "configs/experiment/stage0_flat_deeprvat_full.yaml",
        "matrix_id": MATRIX_ID,
        "source_matrices": sources,
        "merged_matrix": merged,
        "phenotype_table": {
            "path": str(ph_path),
            "n_rows": int(table.num_rows),
            "age_mask_n": int(df["age_mask"].sum()),
            "tissue_mask_n": int(df["tissue_mask"].sum()),
            "sex_mask_n": int(df["sex_mask"].sum()),
        },
        "split": {
            "split_id": split.get("split_id"),
            "mode": split.get("mode"),
            "n_train": len(split.get("train_sample_ids") or []),
            "n_validation": len(split.get("validation_sample_ids") or []),
            "n_external_test": len(split.get("external_test_sample_ids") or []),
            "n_train_studies": len(split.get("train_studies") or []),
            "n_validation_studies": len(split.get("validation_studies") or []),
            "n_external_test_studies": len(split.get("external_test_studies") or []),
        },
        "task": metrics.get("task"),
        "n_classes": metrics.get("n_classes"),
        "n_genes": metrics.get("n_genes"),
        "best_epoch": metrics.get("best_epoch"),
        "best_val_loss": metrics.get("best_val_loss"),
        "final": final,
        "external_test": external,
        "age_standardization": metrics.get("age_standardization"),
        "checkpoints": {
            "best": str(ckpt_dir / "best.pt"),
            "last": str(ckpt_dir / "last.pt"),
            "manifest": str(ckpt_dir / "checkpoint_manifest.json"),
        },
        "tensorboard_dir": str(run_dir / "tb"),
        "metrics_jsonl": str(run_dir / "metrics.jsonl"),
        "pattern": "shared FlatDeepSet + decoupled age/tissue/sex modules + masked loss",
    }
    out_root = project_root / "reports" / "inspection" / "stage0_5d_max_n"
    out_root.mkdir(parents=True, exist_ok=True)
    (out_root / "summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    matrix_lines = [
        f"- `{src['matrix_id']}`: shape `{src['shape']}` "
        f"({src['n_samples']} samples x {src['n_loci']} loci)"
        for src in sources
    ]
    lines = [
        "# Stage 0 Milestone 5d — max-N flat DeepRVAT baseline",
        "",
        "Public model name: **deepMAT**. Shared flat encoder with **decoupled**",
        "age / tissue / sex phenotype modules and **masked** per-trait loss",
        "(DeepRVAT pattern). Disease/cancer not in scope for 5d.",
        "",
        "## Matrices",
        "",
        *matrix_lines,
        f"- **Merged** `{MATRIX_ID}`: shape `{merged['shape']}` "
        f"(GSM-union; notes: {merged.get('notes')})",
        "",
        "## Phenotype table",
        "",
        f"- `{ph_path.name}`: **{payload['phenotype_table']['n_rows']}** samples",
        f"- masks: age={payload['phenotype_table']['age_mask_n']}, "
        f"tissue={payload['phenotype_table']['tissue_mask_n']}, "
        f"sex={payload['phenotype_table']['sex_mask_n']}",
        "",
        "## Train run",
        "",
        f"- run_id: `{RUN_ID}`",
        "- config: `configs/experiment/stage0_flat_deeprvat_full.yaml`",
        f"- split: train={payload['split']['n_train']} / "
        f"val={payload['split']['n_validation']} / "
        f"test={payload['split']['n_external_test']} "
        f"(studies {payload['split']['n_train_studies']} / "
        f"{payload['split']['n_validation_studies']} / "
        f"{payload['split']['n_external_test_studies']})",
        f"- tissue classes: **{metrics.get('n_classes')}**; genes: **{metrics.get('n_genes')}**",
        f"- best_epoch: `{metrics.get('best_epoch')}` "
        f"(best_val_loss={metrics.get('best_val_loss')})",
        f"- external_test tissue accuracy: "
        f"**{external.get('accuracy')}** (n={external.get('tissue_n')})",
        f"- external_test age MAE: **{external.get('mae')}** years (n={external.get('age_n')})",
        f"- external_test sex accuracy: **{external.get('sex_accuracy')}** "
        f"(n={external.get('sex_n')})",
        "",
        "## Artifacts",
        "",
        f"- run: `$MBS_ARTIFACT_ROOT/runs/{RUN_ID}/`",
        f"- checkpoints: `$MBS_ARTIFACT_ROOT/checkpoints/{RUN_ID}/`",
        "- plan: [`docs/plans/milestone-5d-max-n-flat-baseline.md`]"
        "(../../plans/milestone-5d-max-n-flat-baseline.md)",
        "",
    ]
    (out_root / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {out_root}")  # noqa: T201


if __name__ == "__main__":
    main()
