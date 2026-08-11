"""Merge Hub matrix stores into one multitask store (Milestones 5c / 5d)."""

from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from mbs.annotation.manifest import git_commit, sha256_file, utc_now_iso, write_json
from mbs.matrix.store import (
    ARTIFACT_VERSION,
    DEFAULT_DTYPE,
    MATRIX_ORIENTATION,
    MISSING_VALUE_ENCODING,
    matrix_store_paths,
    open_betas_zarr,
    read_locus_index,
    read_sample_index,
    write_betas_zarr,
    write_locus_index,
    write_matrix_manifest,
    write_sample_index,
)
from mbs.training.phenotype_table import (
    build_sample_phenotype_rows,
    build_tissue_ontology,
    default_sex_ontology,
    write_sample_phenotype_table,
    write_sex_ontology,
    write_tissue_ontology,
)


@dataclass(frozen=True, slots=True)
class MultitaskMergeResult:
    matrix_id: str
    output_dir: Path
    n_samples: int
    n_loci: int
    n_deduped: int
    phenotype_table_path: Path
    tissue_ontology_path: Path
    sex_ontology_path: Path | None
    stats: dict[str, Any]


def _assert_identical_loci(*loci: pd.DataFrame) -> None:
    if len(loci) < 2:
        return
    ref = loci[0]
    for other in loci[1:]:
        if len(ref) != len(other):
            raise ValueError(f"locus count mismatch: {len(ref)} vs {len(other)}")
        if not np.array_equal(ref["locus_id"].to_numpy(), other["locus_id"].to_numpy()):
            raise ValueError("locus_id order differs across source matrices; refuse to merge")


def _row_map(sample_index: pd.DataFrame) -> dict[str, int]:
    return {
        str(sid): int(row)
        for sid, row in zip(
            sample_index["sample_id"].astype(str),
            sample_index["row_index"].astype(int),
            strict=True,
        )
    }


def merge_age_tissue_matrices(
    *,
    project_root: Path,
    data_root: Path,
    age_matrix_id: str = "matrix-hub-age-studyholdout-v1",
    tissue_matrix_id: str = "matrix-hub-tissue-studyholdout-v1",
    sex_matrix_id: str | None = None,
    output_matrix_id: str = "matrix-hub-age-tissue-multitask-v1",
    output_dir: Path | None = None,
    phenotype_table_path: Path | None = None,
    tissue_ontology_path: Path | None = None,
    sex_ontology_path: Path | None = None,
    min_tissue_n: int = 10,
) -> MultitaskMergeResult:
    """Stack Hub matrices (identical loci); GSM union; write phenotype table.

    Pack order for beta first-seen: age, tissue, optional sex.
    """
    age_root = data_root / "canonical" / "matrices" / age_matrix_id
    tissue_root = data_root / "canonical" / "matrices" / tissue_matrix_id
    age_paths = matrix_store_paths(age_root)
    tissue_paths = matrix_store_paths(tissue_root)
    required = [
        age_paths.betas_path,
        age_paths.sample_index_path,
        age_paths.locus_index_path,
        tissue_paths.betas_path,
        tissue_paths.sample_index_path,
        tissue_paths.locus_index_path,
    ]
    sex_root: Path | None = None
    sex_paths = None
    if sex_matrix_id is not None:
        sex_root = data_root / "canonical" / "matrices" / sex_matrix_id
        sex_paths = matrix_store_paths(sex_root)
        required.extend(
            [
                sex_paths.betas_path,
                sex_paths.sample_index_path,
                sex_paths.locus_index_path,
            ]
        )
    for path in required:
        if not path.exists():
            raise FileNotFoundError(path)

    age_index = read_sample_index(age_paths.sample_index_path)
    tissue_index = read_sample_index(tissue_paths.sample_index_path)
    age_locus = read_locus_index(age_paths.locus_index_path)
    tissue_locus = read_locus_index(tissue_paths.locus_index_path)
    sex_index = None
    sex_locus = None
    if sex_paths is not None:
        sex_index = read_sample_index(sex_paths.sample_index_path)
        sex_locus = read_locus_index(sex_paths.locus_index_path)
        _assert_identical_loci(age_locus, tissue_locus, sex_locus)
    else:
        _assert_identical_loci(age_locus, tissue_locus)

    age_pheno = pd.read_parquet(age_root / "sample_phenotypes.parquet")
    tissue_pheno = pd.read_parquet(tissue_root / "sample_phenotypes.parquet")
    sex_pheno = (
        pd.read_parquet(sex_root / "sample_phenotypes.parquet") if sex_root is not None else None
    )

    age_betas = open_betas_zarr(age_paths.betas_path)
    tissue_betas = open_betas_zarr(tissue_paths.betas_path)
    sex_betas = open_betas_zarr(sex_paths.betas_path) if sex_paths is not None else None
    n_loci = int(age_locus.shape[0])

    age_rows = _row_map(age_index)
    tissue_rows = _row_map(tissue_index)
    sex_rows = _row_map(sex_index) if sex_index is not None else {}

    age_ids = list(age_rows.keys())
    tissue_only = [sid for sid in tissue_rows if sid not in age_rows]
    sex_only = [sid for sid in sex_rows if sid not in age_rows and sid not in tissue_rows]
    ordered_ids = age_ids + tissue_only + sex_only
    n_sources = len(age_rows) + len(tissue_rows) + len(sex_rows)
    n_deduped = n_sources - len(ordered_ids)

    n_samples = len(ordered_ids)
    betas = np.full((n_samples, n_loci), np.nan, dtype=np.float32)
    for out_i, sid in enumerate(ordered_ids):
        if sid in age_rows:
            betas[out_i, :] = np.asarray(age_betas[age_rows[sid], :], dtype=np.float32)
        elif sid in tissue_rows:
            betas[out_i, :] = np.asarray(tissue_betas[tissue_rows[sid], :], dtype=np.float32)
        else:
            if sex_betas is None:
                raise RuntimeError(f"sample {sid} not in age/tissue and sex betas missing")
            betas[out_i, :] = np.asarray(sex_betas[sex_rows[sid], :], dtype=np.float32)

    out_dir = (
        output_dir
        if output_dir is not None
        else data_root / "canonical" / "matrices" / output_matrix_id
    ).resolve()
    paths = matrix_store_paths(out_dir)
    if paths.root.exists():
        shutil.rmtree(paths.root)
    paths.root.mkdir(parents=True, exist_ok=True)

    write_sample_index(
        paths.sample_index_path,
        sample_ids=ordered_ids,
        source_sample_ids=ordered_ids,
    )
    write_locus_index(
        paths.locus_index_path,
        locus_ids=age_locus["locus_id"].to_numpy(),
        canonical_keys=age_locus["canonical_key"].to_numpy(),
        probe_ids=age_locus["probe_id"].to_numpy(),
    )
    chunks = (min(64, max(1, n_samples)), min(4096, max(1, n_loci)))
    write_betas_zarr(paths.betas_path, betas, chunks=chunks)

    pheno_frames = [("age", age_pheno), ("tissue", tissue_pheno)]
    if sex_pheno is not None:
        pheno_frames.append(("sex", sex_pheno))
    by_frames = {name: frame.set_index("sample_id", drop=False) for name, frame in pheno_frames}
    sidecar_rows: list[dict[str, Any]] = []
    for sid in ordered_ids:
        row = None
        for name, _frame in pheno_frames:
            indexed = by_frames[name]
            if sid in indexed.index:
                row = indexed.loc[sid]
                if isinstance(row, pd.DataFrame):
                    row = row.iloc[0]
                break
        if row is None:
            raise KeyError(f"no phenotype sidecar for sample_id={sid}")
        sidecar_rows.append(row.to_dict())
    pd.DataFrame(sidecar_rows).to_parquet(paths.root / "sample_phenotypes.parquet", index=False)

    source_matrices = [age_matrix_id, tissue_matrix_id]
    if sex_matrix_id is not None:
        source_matrices.append(sex_matrix_id)
    write_json(
        paths.root / "study_subset.json",
        {
            "phenotype_family": "multi",
            "source_matrices": source_matrices,
            "n_samples": n_samples,
            "n_deduped_gsm": n_deduped,
        },
    )

    tissue_labels = []
    for rec in tissue_pheno.to_dict(orient="records"):
        lab = rec.get("phenotype_value") or rec.get("tissue")
        if lab is None or (isinstance(lab, float) and pd.isna(lab)):
            continue
        text = str(lab).strip()
        if text:
            tissue_labels.append(text)
    ontology = build_tissue_ontology(tissue_labels, min_n=min_tissue_n)

    ont_path = (
        tissue_ontology_path
        if tissue_ontology_path is not None
        else data_root / "canonical" / "phenotypes" / "tissue_ontology.yaml"
    )
    write_tissue_ontology(ont_path, ontology)

    sex_ont = default_sex_ontology()
    sex_ont_path: Path | None = None
    if sex_matrix_id is not None:
        sex_ont_path = (
            sex_ontology_path
            if sex_ontology_path is not None
            else data_root / "canonical" / "phenotypes" / "sex_ontology_v1.yaml"
        )
        write_sex_ontology(sex_ont_path, sex_ont)

    sample_index = read_sample_index(paths.sample_index_path)
    table = build_sample_phenotype_rows(
        age_pheno=age_pheno,
        tissue_pheno=tissue_pheno,
        sex_pheno=sex_pheno,
        sample_index=sample_index,
        matrix_id=output_matrix_id,
        ontology=ontology,
        sex_ontology=sex_ont if sex_matrix_id is not None else None,
    )
    default_pheno = (
        "sample_phenotype_table_age_tissue_sex_full_v1.parquet"
        if sex_matrix_id is not None
        else "sample_phenotype_table.parquet"
    )
    pheno_path = (
        phenotype_table_path
        if phenotype_table_path is not None
        else data_root / "canonical" / "phenotypes" / default_pheno
    )
    write_sample_phenotype_table(pheno_path, table)

    sample_index_sha = sha256_file(paths.sample_index_path)
    locus_index_sha = sha256_file(paths.locus_index_path)
    source_files = [
        {
            "path": str(age_paths.manifest_path.resolve()),
            "sha256": sha256_file(age_paths.manifest_path),
            "role": "source_matrix_manifest",
        },
        {
            "path": str(tissue_paths.manifest_path.resolve()),
            "sha256": sha256_file(tissue_paths.manifest_path),
            "role": "source_matrix_manifest",
        },
    ]
    if sex_paths is not None:
        source_files.append(
            {
                "path": str(sex_paths.manifest_path.resolve()),
                "sha256": sha256_file(sex_paths.manifest_path),
                "role": "source_matrix_manifest",
            }
        )
    notes = f"Multitask merge of {' + '.join(source_matrices)}; GSM dedupe={n_deduped}"
    manifest: dict[str, Any] = {
        "artifact_version": ARTIFACT_VERSION,
        "matrix_id": output_matrix_id,
        "study_id": output_matrix_id,
        "platform_id": "HM450",
        "processing_level": "gmqn",
        "genome_build": "GRCh38",
        "shape": [n_samples, n_loci],
        "dtype": DEFAULT_DTYPE,
        "chunks": list(chunks),
        "compression": None,
        "missing_value_encoding": MISSING_VALUE_ENCODING,
        "matrix_orientation": MATRIX_ORIENTATION,
        "matrix_path": str(paths.betas_path),
        "sample_index_path": str(paths.sample_index_path),
        "locus_index_path": str(paths.locus_index_path),
        "sample_index_sha256": sample_index_sha,
        "locus_index_sha256": locus_index_sha,
        "source_files": source_files,
        "conversion_commit": git_commit(project_root),
        "created_at": utc_now_iso(),
        "notes": notes,
    }
    write_matrix_manifest(paths.manifest_path, manifest)

    stats = {
        "n_samples": n_samples,
        "n_loci": n_loci,
        "n_deduped": n_deduped,
        "n_age_mask": int(sum(1 for x in table["age_mask"].tolist() if bool(x))),
        "n_tissue_mask": int(sum(1 for x in table["tissue_mask"].tolist() if bool(x))),
        "n_sex_mask": int(
            sum(1 for x in table["sex_mask"].tolist() if bool(x))
            if "sex_mask" in table.columns
            else 0
        ),
        "tissue_classes": list(ontology.labels),
        "matrix_paths": {
            "root": str(paths.root),
            "betas": str(paths.betas_path),
            "sample_index": str(paths.sample_index_path),
            "locus_index": str(paths.locus_index_path),
            "manifest": str(paths.manifest_path),
        },
    }
    return MultitaskMergeResult(
        matrix_id=output_matrix_id,
        output_dir=paths.root,
        n_samples=n_samples,
        n_loci=n_loci,
        n_deduped=n_deduped,
        phenotype_table_path=pheno_path,
        tissue_ontology_path=ont_path,
        sex_ontology_path=sex_ont_path,
        stats=stats,
    )
