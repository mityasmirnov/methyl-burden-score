"""Convert EWAS Data Hub baseline phenotype packs (wide TSV) to canonical matrices."""

from __future__ import annotations

import hashlib
import zipfile
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, BinaryIO

import numpy as np
import pandas as pd

from mbs.annotation.manifest import git_commit, sha256_file, utc_now_iso, write_json
from mbs.matrix.ewas_db import beta_qc_stats
from mbs.matrix.locus_map import build_probe_locus_map, load_probe_locus_edges
from mbs.matrix.store import (
    ARTIFACT_VERSION,
    DEFAULT_DTYPE,
    MATRIX_ORIENTATION,
    MISSING_VALUE_ENCODING,
    matrix_store_paths,
    write_betas_zarr,
    write_locus_index,
    write_matrix_manifest,
    write_sample_index,
)

# Metadata row keys that appear before probe rows in Hub pack TSVs.
_METADATA_KEYS = frozenset(
    {
        "sample_id",
        "age",
        "tissue",
        "disease",
        "sex",
        "bmi",
        "race",
        "cell_component",
        "sample_type",
        "platform",
        "project_id",
    }
)

_PACK_TXT_NAME = {
    "age": "age_methylation_v1.txt",
    "tissue": "tissue_methylation_v1.txt",
    "disease": "disease_methylation_v1.txt",
    "cancer": "cancer_methylation_v1.txt",
    "blood": "blood_methylation_v1.txt",
    "brain": "brain_methylation_v1.txt",
    "sex": "sex_methylation_v1.txt",
}

_PACK_ZIP_NAME = {
    "age": "age_methylation_v1.zip",
    "tissue": "tissue_methylation_v1.zip",
    "disease": "disease_methylation_v1.zip",
    "cancer": "cancer_methylation_v1.zip",
    "blood": "blood_methylation_v1.zip",
    "brain": "brain_methylation_v1.zip",
    "sex": "sex_methylation_v1.zip",
}


@dataclass(frozen=True, slots=True)
class HubPackConvertResult:
    matrix_id: str
    phenotype_family: str
    output_dir: Path
    sample_ids: tuple[str, ...]
    study_ids: tuple[str, ...]
    stats: dict[str, Any]


def pack_zip_path(data_root: Path, family: str) -> Path:
    if family not in _PACK_ZIP_NAME:
        raise ValueError(f"unsupported Hub pack family: {family}")
    return data_root / "raw" / "ewas_datahub" / "download" / _PACK_ZIP_NAME[family]


def study_ids_from_sample_info(sample_info: pd.DataFrame) -> list[str]:
    """Sorted unique study accessions from sample-info (``study_id`` or ``project_id``)."""
    frame = sample_info
    if "study_id" in frame.columns:
        col = "study_id"
    elif "project_id" in frame.columns:
        col = "project_id"
    else:
        raise ValueError("sample-info requires study_id or project_id")
    values = sorted({str(x) for x in frame[col].tolist() if str(x).strip() and str(x) != "nan"})
    if not values:
        raise ValueError("sample-info has no study accessions")
    return values


def select_samples_for_studies(
    sample_info: pd.DataFrame,
    study_ids: Sequence[str],
    *,
    max_per_study: int | None = None,
) -> pd.DataFrame:
    """Return sample-info rows for the requested studies (stable study then sample order)."""
    wanted = [str(s) for s in study_ids]
    frame = sample_info.copy()
    if "study_id" not in frame.columns and "project_id" in frame.columns:
        frame = frame.copy()
        frame["study_id"] = frame["project_id"]
    frame["study_id"] = frame["study_id"].astype(str)
    frame["sample_id"] = frame["sample_id"].astype(str)
    subset = frame[frame["study_id"].isin(wanted)]
    if subset.empty:
        raise ValueError(f"no sample-info rows for studies={wanted}")
    missing = [s for s in wanted if s not in set(subset["study_id"])]
    if missing:
        raise ValueError(f"studies absent from sample-info: {missing}")
    parts: list[pd.DataFrame] = []
    for study in wanted:
        block = subset[subset["study_id"] == study].sort_values("sample_id")
        if max_per_study is not None:
            block = block.head(int(max_per_study))
        parts.append(block)
    out = pd.concat(parts, ignore_index=True)
    if out["sample_id"].duplicated().any():
        raise ValueError("duplicate sample_id after study selection")
    return out


def _parse_header_sample_ids(header_line: str) -> list[str]:
    parts = header_line.rstrip("\n").split("\t")
    if not parts or parts[0] != "sample_id":
        raise ValueError(f"expected sample_id header row, got {parts[:3]!r}")
    return [str(x) for x in parts[1:]]


def _open_pack_txt(zip_path: Path, family: str) -> tuple[zipfile.ZipFile, BinaryIO, str]:
    zip_path = zip_path.resolve()
    if not zip_path.is_file():
        raise FileNotFoundError(f"Hub profile pack not found: {zip_path}")
    member = _PACK_TXT_NAME[family]
    zf = zipfile.ZipFile(zip_path)
    try:
        fh = zf.open(member)
    except KeyError as exc:
        zf.close()
        raise FileNotFoundError(f"{member} missing inside {zip_path}") from exc
    return zf, fh, member


def stream_pack_betas(
    *,
    zip_path: Path,
    family: str,
    sample_ids: Sequence[str],
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Stream selected columns from a Hub pack TSV into a dense float32 matrix.

    Returns ``(betas [n_samples, n_probes], probe_ids, meta)`` with samples in the
    same order as ``sample_ids``.
    """
    wanted = [str(s) for s in sample_ids]
    zf, fh, member = _open_pack_txt(zip_path, family)
    try:
        header = fh.readline().decode("utf-8", errors="replace")
        pack_samples = _parse_header_sample_ids(header)
        index_by_id = {sid: i for i, sid in enumerate(pack_samples)}
        missing = [sid for sid in wanted if sid not in index_by_id]
        if missing:
            raise KeyError(
                f"{len(missing)} sample_id(s) not in pack {zip_path.name}; first={missing[0]!r}"
            )
        col_idxs = np.asarray([index_by_id[sid] for sid in wanted], dtype=np.int64)
        # Sorted unique pack columns for a single left-to-right scan per probe row.
        order = np.argsort(col_idxs)
        sorted_cols = col_idxs[order]
        inv_order = np.empty_like(order)
        inv_order[order] = np.arange(len(order))

        probe_ids: list[str] = []
        # Accumulate as list of float32 vectors then stack (ponytail: peak RAM ~ probes*samples*4).
        columns: list[np.ndarray] = []
        n_meta_rows = 0
        while True:
            raw = fh.readline()
            if not raw:
                break
            line = raw.decode("utf-8", errors="replace").rstrip("\n")
            if not line:
                continue
            # Split only first field cheaply
            tab = line.find("\t")
            if tab < 0:
                key = line
                rest = ""
            else:
                key = line[:tab]
                rest = line[tab + 1 :]
            if key in _METADATA_KEYS or not key.startswith("cg"):
                n_meta_rows += 1
                continue
            vals_sorted = np.empty(len(wanted), dtype=np.float32)
            start = 0
            col_i = 0
            n_wanted = len(sorted_cols)
            target = int(sorted_cols[col_i])
            field_idx = 0
            while col_i < n_wanted:
                next_tab = rest.find("\t", start)
                if field_idx == target:
                    tok = rest[start:] if next_tab < 0 else rest[start:next_tab]
                    if tok in {"", "NA", "NaN", "nan", "."}:
                        vals_sorted[col_i] = np.nan
                    else:
                        vals_sorted[col_i] = np.float32(tok)
                    col_i += 1
                    if col_i >= n_wanted:
                        break
                    target = int(sorted_cols[col_i])
                if next_tab < 0:
                    break
                start = next_tab + 1
                field_idx += 1
            if col_i != n_wanted:
                raise ValueError(
                    f"probe {key}: could not extract selected columns "
                    f"(got {col_i}/{n_wanted}; pack_width={len(pack_samples)})"
                )
            vals = vals_sorted[inv_order]
            probe_ids.append(key)
            columns.append(vals)

        if not probe_ids:
            raise ValueError(f"no probe rows found in {member}")
        # columns are per-probe length n_samples → stack to [n_probes, n_samples] then T
        probe_mat = np.stack(columns, axis=0)
        betas = np.ascontiguousarray(probe_mat.T)
        meta = {
            "pack_member": member,
            "n_pack_samples": len(pack_samples),
            "n_meta_rows_skipped": n_meta_rows,
            "n_probes": len(probe_ids),
            "n_selected_samples": len(wanted),
        }
        return betas, np.asarray(probe_ids, dtype=object), meta
    finally:
        fh.close()
        zf.close()


def convert_hub_pack_subset(
    *,
    project_root: Path,
    data_root: Path,
    annotations_dir: Path,
    phenotype_family: str,
    study_ids: Sequence[str],
    matrix_id: str,
    output_dir: Path,
    platform_id: str = "HM450",
    processing_level: str = "gmqn",
    max_per_study: int | None = None,
    sample_info_path: Path | None = None,
) -> HubPackConvertResult:
    """Convert a study-subset of one Hub baseline pack into a canonical matrix store."""
    family = phenotype_family
    if family not in _PACK_ZIP_NAME:
        raise ValueError(f"unsupported phenotype_family: {family}")
    zip_path = pack_zip_path(data_root, family)
    info_path = (
        sample_info_path
        if sample_info_path is not None
        else data_root / "canonical" / "phenotypes" / f"{family}_sample_info.parquet"
    )
    if not info_path.is_file():
        raise FileNotFoundError(f"sample-info parquet required before pack convert: {info_path}")
    sample_info = pd.read_parquet(info_path)
    selected = select_samples_for_studies(sample_info, study_ids, max_per_study=max_per_study)
    sample_ids = selected["sample_id"].astype(str).tolist()
    study_list = [str(s) for s in study_ids]

    raw_betas, probe_ids, pack_meta = stream_pack_betas(
        zip_path=zip_path, family=family, sample_ids=sample_ids
    )
    # raw_betas is [n_samples, n_probes]; map probes → loci
    edges = load_probe_locus_edges(annotations_dir, platform_id=platform_id)
    locus_map = build_probe_locus_map(probe_ids, edges, platform_id=platform_id)
    probe_to_col = {str(pid): i for i, pid in enumerate(locus_map.probe_ids)}

    n_samples = len(sample_ids)
    n_loci = len(locus_map.locus_ids)
    betas = np.full((n_samples, n_loci), np.nan, dtype=np.float32)
    for p_idx, pid in enumerate(probe_ids.astype(str)):
        col = probe_to_col.get(pid)
        if col is None:
            continue
        betas[:, col] = raw_betas[:, p_idx]

    del raw_betas

    n_out_of_range = 0
    n_missing_cells = 0
    finite_chunks: list[np.ndarray] = []
    for row in range(n_samples):
        qc = beta_qc_stats(betas[row, :])
        n_out_of_range += int(qc["n_out_of_range"])
        n_missing_cells += int(qc["n_missing"])
        mask = np.isfinite(betas[row, :])
        if mask.any():
            finite_chunks.append(betas[row, mask].copy())

    output_dir = output_dir.resolve()
    paths = matrix_store_paths(output_dir)
    paths.root.mkdir(parents=True, exist_ok=True)
    write_sample_index(
        paths.sample_index_path,
        sample_ids=sample_ids,
        source_sample_ids=sample_ids,
    )
    # Sidecar study membership for study-grouped splits
    keep_cols = [
        c
        for c in (
            "sample_id",
            "study_id",
            "platform",
            "phenotype_value",
            "phenotype_value_numeric",
            "sample_type",
            "tissue",
            "disease",
            "age",
            "sex",
        )
        if c in selected.columns
    ]
    study_frame = selected[keep_cols].copy()
    study_frame.to_parquet(paths.root / "sample_phenotypes.parquet", index=False)
    write_json(
        paths.root / "study_subset.json",
        {
            "phenotype_family": family,
            "study_ids": study_list,
            "n_samples": n_samples,
            "max_per_study": max_per_study,
            "source_pack": str(zip_path),
        },
    )
    write_locus_index(
        paths.locus_index_path,
        locus_ids=locus_map.locus_ids,
        canonical_keys=locus_map.canonical_keys,
        probe_ids=locus_map.probe_ids,
    )
    chunks = (min(64, max(1, n_samples)), min(4096, max(1, n_loci)))
    write_betas_zarr(paths.betas_path, betas, chunks=chunks)

    if finite_chunks:
        all_finite = np.concatenate(finite_chunks)
        beta_min = float(all_finite.min())
        beta_max = float(all_finite.max())
        beta_mean = float(all_finite.mean())
        n_finite = int(all_finite.size)
    else:
        beta_min = float("nan")
        beta_max = float("nan")
        beta_mean = float("nan")
        n_finite = 0

    sample_index_sha = sha256_file(paths.sample_index_path)
    locus_index_sha = sha256_file(paths.locus_index_path)
    notes = (
        f"EWAS Data Hub baseline pack {family}; studies={','.join(study_list)}; "
        f"platform={platform_id}; unmapped_probes={len(locus_map.unmapped_probe_ids)}; "
        f"collapsed_probes={locus_map.n_collapsed_probes}; pack={zip_path.name}"
    )
    # Manifest study_id: multi-study packs use matrix_id as the cohort token.
    if zip_path.stat().st_size >= 50_000_000:
        digest = hashlib.sha256(
            f"{zip_path.name}:{zip_path.stat().st_size}:{pack_meta['pack_member']}".encode()
        ).hexdigest()
        source_files: list[dict[str, Any]] = [
            {
                "path": str(zip_path),
                "sha256": digest,
                "role": "hub_baseline_pack",
                "bytes": zip_path.stat().st_size,
                "sha256_note": "content-address of name+size (pack too large to hash)",
            }
        ]
    else:
        source_files = [
            {
                "path": str(zip_path),
                "sha256": sha256_file(zip_path),
                "role": "hub_baseline_pack",
            }
        ]

    manifest: dict[str, Any] = {
        "artifact_version": ARTIFACT_VERSION,
        "matrix_id": matrix_id,
        "study_id": matrix_id,
        "platform_id": platform_id,
        "processing_level": processing_level,
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
        "n_study_loci": n_loci,
        "n_observed_probes": locus_map.n_observed_probes,
        "n_mapped_probes": locus_map.n_mapped_probes,
        "n_unmapped_probes": len(locus_map.unmapped_probe_ids),
        "n_collapsed_probes": locus_map.n_collapsed_probes,
        "n_finite_betas": n_finite,
        "n_missing_cells": n_missing_cells,
        "n_out_of_range": n_out_of_range,
        "beta_min": beta_min,
        "beta_max": beta_max,
        "beta_mean": beta_mean,
        "study_ids": study_list,
        "pack_meta": pack_meta,
        "matrix_paths": {
            "root": str(paths.root),
            "betas": str(paths.betas_path),
            "sample_index": str(paths.sample_index_path),
            "locus_index": str(paths.locus_index_path),
            "manifest": str(paths.manifest_path),
            "sample_phenotypes": str(paths.root / "sample_phenotypes.parquet"),
        },
    }
    write_json(paths.root / "conversion_stats.json", stats)
    return HubPackConvertResult(
        matrix_id=matrix_id,
        phenotype_family=family,
        output_dir=paths.root,
        sample_ids=tuple(sample_ids),
        study_ids=tuple(study_list),
        stats=stats,
    )
