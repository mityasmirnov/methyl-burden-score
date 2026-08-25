"""Convert EWAS Data Hub baseline phenotype packs (wide TSV) to canonical matrices."""

from __future__ import annotations

import zipfile
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, BinaryIO

import numpy as np
import pandas as pd

from mbs.annotation.manifest import git_commit, sha256_file, utc_now_iso, write_json
from mbs.matrix.ewas_db import beta_qc_stats
from mbs.matrix.locus_map import (
    COLLAPSE_IDENTITY,
    COLLAPSE_MEAN,
    COLLAPSE_MEDIAN,
    ProbeLocusMap,
    build_probe_locus_map,
    load_probe_locus_edges,
)
from mbs.matrix.store import (
    ARTIFACT_VERSION,
    DEFAULT_DTYPE,
    DEFAULT_ZARR_COMPRESSION,
    MATRIX_ORIENTATION,
    MISSING_VALUE_ENCODING,
    create_betas_zarr,
    matrix_store_paths,
    open_betas_zarr,
    write_locus_index,
    write_matrix_manifest,
    write_sample_index,
)
from mbs.platform_id import normalize_platform

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
    "bmi": "bmi_methylation_v1.txt",
    "ancestry": "ancestry_category_methylation_v1.txt",
}

_PACK_ZIP_NAME = {
    "age": "age_methylation_v1.zip",
    "tissue": "tissue_methylation_v1.zip",
    "disease": "disease_methylation_v1.zip",
    "cancer": "cancer_methylation_v1.zip",
    "blood": "blood_methylation_v1.zip",
    "brain": "brain_methylation_v1.zip",
    "sex": "sex_methylation_v1.zip",
    "bmi": "bmi_methylation_v1.zip",
    "ancestry": "ancestry_category_methylation_v1.zip",
}

SUPPORTED_PACK_FAMILIES: tuple[str, ...] = tuple(_PACK_ZIP_NAME.keys())

_PHENOTYPE_KEEP_COLS = (
    "sample_id",
    "study_id",
    "platform",
    "phenotype_value",
    "phenotype_value_numeric",
    "phenotype_family",
    "sample_type",
    "tissue",
    "disease",
    "age",
    "sex",
    "bmi",
    "race",
    "cell_component",
)


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


def pack_txt_name(family: str) -> str:
    if family not in _PACK_TXT_NAME:
        raise ValueError(f"unsupported Hub pack family: {family}")
    return _PACK_TXT_NAME[family]


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
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return ``(unique_gsm_frame, long_form_frame)`` for requested studies.

    Matrix rows use unique ``sample_id`` (stable study then sample order).
    ``long_form_frame`` keeps every sample-info row (disease/cancer multi-label).
    """
    wanted = [str(s) for s in study_ids]
    frame = sample_info.copy()
    if "study_id" not in frame.columns and "project_id" in frame.columns:
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
            # Cap unique GSMs per study, keep all long-form rows for those GSMs.
            unique_ids = list(dict.fromkeys(block["sample_id"].tolist()))
            keep_ids = set(unique_ids[: int(max_per_study)])
            block = block[block["sample_id"].isin(keep_ids)]
        parts.append(block)
    long_form = pd.concat(parts, ignore_index=True)

    if "platform" in long_form.columns:
        plat = long_form.assign(platform=long_form["platform"].astype(str))
        n_plat = plat.groupby("sample_id")["platform"].nunique()
        conflicts = n_plat[n_plat > 1]
        if not conflicts.empty:
            sid = str(conflicts.index[0])
            vals = sorted(plat.loc[plat["sample_id"] == sid, "platform"].unique())
            raise ValueError(f"conflicting platform for sample_id={sid!r}: {vals}")

    # Unique GSM order: first appearance in long_form (stable study then sample).
    unique = long_form.drop_duplicates(subset=["sample_id"], keep="first").reset_index(drop=True)
    return unique, long_form


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


_MISSING_TOKS = frozenset({b"", b"NA", b"NaN", b"nan", b"."})


def _extract_selected_values(
    rest: bytes,
    *,
    sorted_cols: np.ndarray,
    inv_order: np.ndarray,
) -> np.ndarray:
    """Parse selected sample columns from one TSV probe value line (bytes)."""
    if rest.endswith(b"\n"):
        rest = rest[:-1]
    if rest.endswith(b"\r"):
        rest = rest[:-1]
    n_wanted = len(sorted_cols)
    vals_sorted = np.empty(n_wanted, dtype=np.float32)
    start = 0
    col_i = 0
    target = int(sorted_cols[col_i])
    field_idx = 0
    while col_i < n_wanted:
        next_tab = rest.find(b"\t", start)
        if field_idx == target:
            tok = rest[start:] if next_tab < 0 else rest[start:next_tab]
            if tok in _MISSING_TOKS:
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
        raise ValueError(f"could not extract selected columns (got {col_i}/{n_wanted})")
    return vals_sorted[inv_order]


def scan_pack_probe_ids(
    *,
    zip_path: Path,
    family: str,
    sample_ids: Sequence[str],
) -> tuple[np.ndarray, dict[str, Any], np.ndarray, np.ndarray]:
    """Pass 1: collect probe IDs and column reorder maps (no beta storage)."""
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
        order = np.argsort(col_idxs)
        sorted_cols = col_idxs[order]
        inv_order = np.empty_like(order)
        inv_order[order] = np.arange(len(order))

        probe_ids: list[str] = []
        n_meta_rows = 0
        while True:
            raw = fh.readline()
            if not raw:
                break
            if raw in {b"\n", b"\r\n"}:
                continue
            tab = raw.find(b"\t")
            key_b = raw.rstrip(b"\r\n") if tab < 0 else raw[:tab]
            key = key_b.decode("utf-8", errors="replace")
            if key in _METADATA_KEYS or not key.startswith("cg"):
                n_meta_rows += 1
                continue
            probe_ids.append(key)

        if not probe_ids:
            raise ValueError(f"no probe rows found in {member}")
        meta = {
            "pack_member": member,
            "n_pack_samples": len(pack_samples),
            "n_meta_rows_skipped": n_meta_rows,
            "n_probes": len(probe_ids),
            "n_selected_samples": len(wanted),
        }
        return (
            np.asarray(probe_ids, dtype=object),
            meta,
            sorted_cols,
            inv_order,
        )
    finally:
        fh.close()
        zf.close()


def stream_pack_betas(
    *,
    zip_path: Path,
    family: str,
    sample_ids: Sequence[str],
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Stream selected columns into a dense float32 matrix (small packs / tests).

    Prefer ``convert_hub_pack_subset`` for full packs (streams directly to Zarr).
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
        order = np.argsort(col_idxs)
        sorted_cols = col_idxs[order]
        inv_order = np.empty_like(order)
        inv_order[order] = np.arange(len(order))

        probe_ids: list[str] = []
        columns: list[np.ndarray] = []
        n_meta_rows = 0
        while True:
            raw = fh.readline()
            if not raw:
                break
            if raw in {b"\n", b"\r\n"}:
                continue
            tab = raw.find(b"\t")
            key_b = raw.rstrip(b"\r\n") if tab < 0 else raw[:tab]
            key = key_b.decode("ascii")
            if key in _METADATA_KEYS or not key.startswith("cg"):
                n_meta_rows += 1
                continue
            rest = b"" if tab < 0 else raw[tab + 1 :]
            try:
                vals = _extract_selected_values(rest, sorted_cols=sorted_cols, inv_order=inv_order)
            except ValueError as exc:
                raise ValueError(f"probe {key}: {exc}; pack_width={len(pack_samples)}") from exc
            probe_ids.append(key)
            columns.append(vals)

        if not probe_ids:
            raise ValueError(f"no probe rows found in {member}")
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


def _probe_to_columns(locus_map: ProbeLocusMap) -> dict[str, list[int]]:
    """Map each contributing probe_id → list of matrix column indices."""
    out: dict[str, list[int]] = {}
    for col, probes in enumerate(locus_map.contributing_probe_ids):
        for pid in probes:
            out.setdefault(pid, []).append(col)
    return out


def _collapse_buffered(
    values: list[np.ndarray],
    method: str,
) -> np.ndarray:
    stacked = np.stack(values, axis=0)
    if method == COLLAPSE_MEDIAN:
        return np.nanmedian(stacked, axis=0).astype(np.float32, copy=False)
    if method == COLLAPSE_MEAN:
        with np.errstate(all="ignore"):
            return np.nanmean(stacked, axis=0).astype(np.float32, copy=False)
    return stacked[0].astype(np.float32, copy=False)


def _stream_pack_to_zarr(
    *,
    zip_path: Path,
    family: str,
    sample_ids: Sequence[str],
    locus_map: ProbeLocusMap,
    betas_path: Path,
    chunks: tuple[int, int],
) -> dict[str, Any]:
    """Pass 2: stream probe rows into a scratch memmap, then compressed Zarr.

    Peak RAM is one probe vector plus small multi-probe buffers — not a dense
    ``[n_samples, n_probes]`` stack. Scratch memmap lives next to the output.
    """
    wanted = [str(s) for s in sample_ids]
    n_samples = len(wanted)
    n_loci = len(locus_map.locus_ids)
    probe_cols = _probe_to_columns(locus_map)
    multi_cols = {
        i for i, method in enumerate(locus_map.collapse_method) if method != COLLAPSE_IDENTITY
    }

    scratch = betas_path.parent / ".betas_scratch.f32"
    if scratch.exists():
        scratch.unlink()
    # Locus-major memmap so each probe write is contiguous; transpose into Zarr.
    mm = np.memmap(scratch, dtype=np.float32, mode="w+", shape=(n_loci, n_samples))
    written = np.zeros(n_loci, dtype=bool)

    # Multi-probe buffers only (few HM450 collisions) — not a full dense stack.
    multi_buf: dict[int, list[np.ndarray]] = {c: [] for c in multi_cols}
    mean_sum: dict[int, np.ndarray] = {}
    mean_count: dict[int, np.ndarray] = {}

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
        order = np.argsort(col_idxs)
        sorted_cols = col_idxs[order]
        inv_order = np.empty_like(order)
        inv_order[order] = np.arange(len(order))

        n_meta_rows = 0
        n_probes_seen = 0
        while True:
            raw = fh.readline()
            if not raw:
                break
            if raw in {b"\n", b"\r\n"}:
                continue
            tab = raw.find(b"\t")
            key_b = raw.rstrip(b"\r\n") if tab < 0 else raw[:tab]
            key = key_b.decode("ascii")
            if key in _METADATA_KEYS or not key.startswith("cg"):
                n_meta_rows += 1
                continue
            cols = probe_cols.get(key)
            if cols is None:
                continue
            rest = b"" if tab < 0 else raw[tab + 1 :]
            try:
                vals = _extract_selected_values(rest, sorted_cols=sorted_cols, inv_order=inv_order)
            except ValueError as exc:
                raise ValueError(f"probe {key}: {exc}; pack_width={len(pack_samples)}") from exc
            n_probes_seen += 1
            for col in cols:
                method = locus_map.collapse_method[col]
                if method == COLLAPSE_IDENTITY:
                    mm[col, :] = vals
                    written[col] = True
                elif method == COLLAPSE_MEAN:
                    if col not in mean_sum:
                        mean_sum[col] = np.zeros(n_samples, dtype=np.float64)
                        mean_count[col] = np.zeros(n_samples, dtype=np.int32)
                    finite = np.isfinite(vals)
                    mean_sum[col][finite] += vals[finite]
                    mean_count[col][finite] += 1
                else:
                    multi_buf[col].append(vals.copy())
    finally:
        fh.close()
        zf.close()

    for col, total in mean_sum.items():
        count = mean_count[col]
        out = np.full(n_samples, np.nan, dtype=np.float32)
        ok = count > 0
        out[ok] = (total[ok] / count[ok]).astype(np.float32)
        mm[col, :] = out
        written[col] = True

    for col, buf in multi_buf.items():
        if not buf:
            continue
        method = locus_map.collapse_method[col]
        mm[col, :] = _collapse_buffered(buf, method)
        written[col] = True

    missing_cols = np.flatnonzero(~written)
    if missing_cols.size:
        mm[missing_cols, :] = np.nan
    mm.flush()

    array = create_betas_zarr(
        betas_path,
        n_samples=n_samples,
        n_loci=n_loci,
        chunks=chunks,
    )
    row_chunk = chunks[0]
    for start in range(0, n_samples, row_chunk):
        stop = min(start + row_chunk, n_samples)
        array[start:stop, :] = np.asarray(mm[:, start:stop]).T
    del mm
    scratch.unlink(missing_ok=True)

    return {
        "pack_member": member,
        "n_pack_samples": len(pack_samples),
        "n_meta_rows_skipped": n_meta_rows,
        "n_probes": n_probes_seen,
        "n_selected_samples": n_samples,
        "n_multi_probe_loci": len(multi_cols),
    }


def _manifest_platform_id(unique_samples: pd.DataFrame, fallback: str) -> str:
    if "platform" not in unique_samples.columns:
        return normalize_platform(fallback) or fallback
    values = sorted(
        {
            normalize_platform(p) or str(p).strip()
            for p in unique_samples["platform"].tolist()
            if str(p).strip() and str(p) != "nan"
        }
    )
    if not values:
        return normalize_platform(fallback) or fallback
    if len(values) == 1:
        return values[0]
    return "mixed"


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
    unique, long_form = select_samples_for_studies(
        sample_info, study_ids, max_per_study=max_per_study
    )
    sample_ids = unique["sample_id"].astype(str).tolist()
    study_list = [str(s) for s in study_ids]

    # Pass 1: probe IDs only (no dense beta stack).
    probe_ids, scan_meta, _sorted_cols, _inv = scan_pack_probe_ids(
        zip_path=zip_path, family=family, sample_ids=sample_ids
    )
    edges = load_probe_locus_edges(annotations_dir, platform_id=platform_id)
    locus_map = build_probe_locus_map(probe_ids, edges, platform_id=platform_id)

    n_samples = len(sample_ids)
    n_loci = len(locus_map.locus_ids)
    chunks = (min(64, max(1, n_samples)), min(4096, max(1, n_loci)))

    output_dir = output_dir.resolve()
    paths = matrix_store_paths(output_dir)
    paths.root.mkdir(parents=True, exist_ok=True)

    pack_meta = _stream_pack_to_zarr(
        zip_path=zip_path,
        family=family,
        sample_ids=sample_ids,
        locus_map=locus_map,
        betas_path=paths.betas_path,
        chunks=chunks,
    )
    pack_meta = {**scan_meta, **pack_meta}

    write_sample_index(
        paths.sample_index_path,
        sample_ids=sample_ids,
        source_sample_ids=sample_ids,
    )
    keep_cols = [c for c in _PHENOTYPE_KEEP_COLS if c in long_form.columns]
    study_frame = long_form[keep_cols].copy()
    study_frame.to_parquet(paths.root / "sample_phenotypes.parquet", index=False)
    write_json(
        paths.root / "study_subset.json",
        {
            "phenotype_family": family,
            "study_ids": study_list,
            "n_samples": n_samples,
            "n_phenotype_rows": len(long_form),
            "max_per_study": max_per_study,
            "source_pack": str(zip_path),
        },
    )
    write_locus_index(
        paths.locus_index_path,
        locus_ids=locus_map.locus_ids,
        canonical_keys=locus_map.canonical_keys,
        probe_ids=locus_map.probe_ids,
        annotation_status=locus_map.annotation_status,
        contributing_probe_ids=locus_map.contributing_probe_ids,
        collapse_method=locus_map.collapse_method,
    )

    # QC over sample chunks (avoid loading full matrix when huge).
    betas = open_betas_zarr(paths.betas_path)
    n_out_of_range = 0
    n_missing_cells = 0
    finite_chunks: list[np.ndarray] = []
    row_chunk = min(64, max(1, n_samples))
    for start in range(0, n_samples, row_chunk):
        stop = min(start + row_chunk, n_samples)
        block = np.asarray(betas[start:stop, :], dtype=np.float32)
        for row in range(block.shape[0]):
            qc = beta_qc_stats(block[row, :])
            n_out_of_range += int(qc["n_out_of_range"])
            n_missing_cells += int(qc["n_missing"])
            mask = np.isfinite(block[row, :])
            if mask.any():
                # Subsample finite values for global min/max/mean (ponytail: full concat is huge).
                finite_chunks.append(block[row, mask][:: max(1, mask.sum() // 10_000)].copy())

    if finite_chunks:
        all_finite = np.concatenate(finite_chunks)
        beta_min = float(all_finite.min())
        beta_max = float(all_finite.max())
        beta_mean = float(all_finite.mean())
        n_finite = int(n_samples * n_loci - n_missing_cells)
    else:
        beta_min = float("nan")
        beta_max = float("nan")
        beta_mean = float("nan")
        n_finite = 0

    sample_index_sha = sha256_file(paths.sample_index_path)
    locus_index_sha = sha256_file(paths.locus_index_path)
    manifest_platform = _manifest_platform_id(unique, platform_id)
    notes = (
        f"EWAS Data Hub baseline pack {family}; studies={','.join(study_list)}; "
        f"probe_edges_platform={platform_id}; manifest_platform={manifest_platform}; "
        f"unmapped_probes={len(locus_map.unmapped_probe_ids)}; "
        f"residual_probes={locus_map.n_residual_probes}; "
        f"collapsed_probes={locus_map.n_collapsed_probes}; pack={zip_path.name}"
    )
    digest = sha256_file(zip_path)
    source_files: list[dict[str, Any]] = [
        {
            "path": str(zip_path.resolve()),
            "sha256": digest,
            "byte_size": int(zip_path.stat().st_size),
        }
    ]

    manifest: dict[str, Any] = {
        "artifact_version": ARTIFACT_VERSION,
        "matrix_id": matrix_id,
        "study_id": matrix_id,
        "platform_id": manifest_platform,
        "processing_level": processing_level,
        "genome_build": "GRCh38",
        "shape": [n_samples, n_loci],
        "dtype": DEFAULT_DTYPE,
        "chunks": list(chunks),
        "compression": DEFAULT_ZARR_COMPRESSION,
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
        "n_phenotype_rows": len(long_form),
        "n_study_loci": n_loci,
        "n_observed_probes": locus_map.n_observed_probes,
        "n_mapped_probes": locus_map.n_mapped_probes,
        "n_unmapped_probes": len(locus_map.unmapped_probe_ids),
        "n_residual_probes": locus_map.n_residual_probes,
        "n_collapsed_probes": locus_map.n_collapsed_probes,
        "n_finite_betas": n_finite,
        "n_missing_cells": n_missing_cells,
        "n_out_of_range": n_out_of_range,
        "beta_min": beta_min,
        "beta_max": beta_max,
        "beta_mean": beta_mean,
        "study_ids": study_list,
        "pack_meta": pack_meta,
        "probe_edges_platform": platform_id,
        "manifest_platform": manifest_platform,
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
