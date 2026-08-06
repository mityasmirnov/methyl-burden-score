"""Export CpGPT2M sequence-adapter embeddings for the canonical locus registry."""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from importlib import import_module
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import torch

from mbs.annotation.manifest import git_commit, sha256_file, write_json
from mbs.static_features.coordinates import mbs_locus_to_cpgpt_location
from mbs.static_features.cpgpt_adapter import load_small_sequence_adapter
from mbs.static_features.manifest import (
    ARTIFACT_VERSION,
    SOURCE_REPOSITORY,
    utc_now_iso,
)
from mbs.static_features.store import (
    static_feature_store_paths,
    write_artifact,
    write_embeddings_zarr,
    write_loci_index,
)
from mbs.static_features.validate_export import (
    embedding_summary_stats,
    validate_embeddings_array,
    validate_loci_frame,
)

DEFAULT_FEATURE_SET_ID = "cpgpt2m_adapter_128_v1"
DNA_LLM = "nucleotide-transformer-v2-500m-multi-species"
SPECIES = "homo_sapiens"
CONTEXT_LENGTH = 2001
INPUT_DIMENSION = 1024
OUTPUT_DIMENSION = 128
STORAGE_DTYPE = "float16"
DEFAULT_BATCH_SIZE = 8192


@dataclass(frozen=True, slots=True)
class CpGPTExportResult:
    feature_set_id: str
    output_dir: Path
    report_dir: Path | None
    manifest: dict[str, Any]
    stats: dict[str, Any]


def _require_cpgpt_downloads() -> None:
    try:
        import_module("cpgpt.downloads")
    except ImportError as error:
        raise RuntimeError(
            "CpGPT is not installed. Install the optional extra first:\n"
            "  uv sync --all-groups --extra cpgpt"
        ) from error


def _load_location_index(dependencies_path: Path) -> dict[str, int]:
    # sqlitedict ships with the optional cpgpt extra, not the default MBS env.
    from sqlitedict import SqliteDict  # noqa: PLC0415

    metadata_path = dependencies_path / "ensembl_metadata.db"
    if not metadata_path.is_file():
        raise FileNotFoundError(f"CpGPT ensembl metadata not found: {metadata_path}")
    with SqliteDict(str(metadata_path), autocommit=False) as db:
        if SPECIES not in db:
            raise KeyError(f"{SPECIES} missing from {metadata_path}")
        species_meta = db[SPECIES]
    if DNA_LLM not in species_meta:
        raise KeyError(f"{DNA_LLM} missing from ensembl metadata for {SPECIES}")
    if CONTEXT_LENGTH not in species_meta[DNA_LLM]:
        raise KeyError(
            f"context_length={CONTEXT_LENGTH} missing from ensembl metadata for {DNA_LLM}"
        )
    loc_index = species_meta[DNA_LLM][CONTEXT_LENGTH]
    if not isinstance(loc_index, dict) or not loc_index:
        raise ValueError("CpGPT location index is empty or invalid")
    return {str(key): int(value) for key, value in loc_index.items()}


def _open_dna_mmap(dependencies_path: Path, *, max_row_index: int) -> np.memmap:
    mmap_path = (
        dependencies_path
        / "dna_embeddings"
        / SPECIES
        / DNA_LLM
        / f"{CONTEXT_LENGTH}bp_dna_embeddings.mmap"
    )
    if not mmap_path.exists():
        raise FileNotFoundError(f"CpGPT DNA embedding mmap not found: {mmap_path}")
    file_rows = mmap_path.stat().st_size // (INPUT_DIMENSION * 4)
    if file_rows < 1:
        raise ValueError(f"DNA embedding mmap is empty: {mmap_path}")
    if max_row_index >= file_rows:
        raise IndexError(f"source embedding row {max_row_index} exceeds mmap rows {file_rows}")
    return np.memmap(
        mmap_path,
        dtype=np.float32,
        mode="r",
        shape=(file_rows, INPUT_DIMENSION),
    )


def _encode_batches(
    adapter: Any,
    dna_rows: np.ndarray,
    *,
    device: str,
    batch_size: int,
) -> np.ndarray:
    n = int(dna_rows.shape[0])
    out = np.empty((n, OUTPUT_DIMENSION), dtype=np.float16)
    torch_device = torch.device(device)
    with torch.inference_mode():
        for start in range(0, n, batch_size):
            end = min(start + batch_size, n)
            batch = torch.as_tensor(
                dna_rows[start:end],
                dtype=torch.float32,
                device=torch_device,
            ).unsqueeze(1)
            encoded = adapter(batch).squeeze(1)
            out[start:end] = encoded.detach().cpu().numpy().astype(np.float16, copy=False)
    return out


def _build_loci_frame(
    loci: pd.DataFrame,
    *,
    location_index: dict[str, int],
) -> tuple[pd.DataFrame, np.ndarray]:
    source_keys: list[str | None] = []
    source_rows: list[int | None] = []
    statuses: list[str] = []
    mapped_source_rows: list[int] = []

    for chrom, pos in zip(loci["chromosome"], loci["position"], strict=True):
        key = mbs_locus_to_cpgpt_location(str(chrom), int(pos))
        if key is None:
            source_keys.append(None)
            source_rows.append(None)
            statuses.append("missing")
            continue
        row = location_index.get(key)
        if row is None:
            source_keys.append(key)
            source_rows.append(None)
            statuses.append("missing")
            continue
        source_keys.append(key)
        source_rows.append(row)
        statuses.append("mapped")
        mapped_source_rows.append(row)

    embedding_row: list[int | None] = [None] * len(loci)
    mapped_i = 0
    for i, status in enumerate(statuses):
        if status == "mapped":
            embedding_row[i] = mapped_i
            mapped_i += 1

    frame = pd.DataFrame(
        {
            "embedding_row": pd.array(embedding_row, dtype="Int64"),
            "locus_id": loci["locus_id"].astype("uint64"),
            "canonical_key": loci["canonical_key"].astype("string"),
            "source_location_key": pd.array(source_keys, dtype="string"),
            "source_embedding_row": pd.array(source_rows, dtype="Int64"),
            "mapping_status": pd.array(statuses, dtype="string"),
        }
    )
    return frame, np.asarray(mapped_source_rows, dtype=np.int64)


def _write_report(
    report_dir: Path,
    *,
    manifest: dict[str, Any],
    stats: dict[str, Any],
) -> None:
    report_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "feature_set_id": manifest["feature_set_id"],
        "manifest": manifest,
        "stats": stats,
    }
    write_json(report_dir / "summary.json", payload)
    lines = [
        "# CpGPT static feature export",
        "",
        f"- feature_set_id: `{manifest['feature_set_id']}`",
        f"- genome_build: `{manifest['genome_build']}`",
        f"- source_commit: `{manifest['source_commit']}`",
        f"- checkpoint_sha256: `{manifest['checkpoint_sha256']}`",
        f"- locus_table_sha256: `{manifest['locus_table_sha256']}`",
        f"- output_dimension: `{manifest['output_dimension']}`",
        f"- storage_dtype: `{manifest['storage_dtype']}`",
        f"- n_loci (registry): `{stats['n_loci']}`",
        f"- n_mapped: `{stats['n_mapped']}`",
        f"- n_missing: `{stats['n_missing']}`",
        f"- mapping_rate: `{stats['mapping_rate']:.6f}`",
        f"- norm_mean: `{stats['embedding_stats']['norm_mean']:.6f}`",
        f"- dim_var_mean: `{stats['embedding_stats']['dim_var_mean']:.6f}`",
        f"- export_command: `{manifest['export_command']}`",
        "",
    ]
    (report_dir / "summary.md").write_text("\n".join(lines), encoding="utf-8")


def export_cpgpt_adapter(
    *,
    project_root: Path,
    loci_path: Path,
    annotations_manifest_path: Path | None = None,
    output_dir: Path,
    report_dir: Path | None = None,
    feature_set_id: str = DEFAULT_FEATURE_SET_ID,
    device: str = "cuda",
    batch_size: int = DEFAULT_BATCH_SIZE,
    export_command: str | None = None,
    cache_dir: Path | None = None,
) -> CpGPTExportResult:
    """Export CpGPT sequence-adapter embeddings for every registry locus."""
    _require_cpgpt_downloads()
    from cpgpt.downloads import download_cpgpt  # noqa: PLC0415

    resolved_loci = loci_path.resolve()
    if not resolved_loci.is_file():
        raise FileNotFoundError(f"loci parquet not found: {resolved_loci}")

    loci_sha = sha256_file(resolved_loci)
    if annotations_manifest_path is not None:
        ann_manifest = json.loads(annotations_manifest_path.read_text(encoding="utf-8"))
        expected = ann_manifest.get("loci_sha256")
        if expected is None:
            raise ValueError(
                f"annotations manifest missing loci_sha256: {annotations_manifest_path}"
            )
        if expected != loci_sha:
            raise ValueError(f"locus table hash mismatch: file={loci_sha} manifest={expected}")

    hf_home = Path(cache_dir or os.environ.get("HF_HOME", project_root / "cache" / "huggingface"))
    resources = download_cpgpt(model="small", species="human", cache_dir=str(hf_home))
    if resources.checkpoint_path is None or resources.config_path is None:
        raise RuntimeError("CpGPT small checkpoint/config not resolved from cache")
    if resources.dependencies_path is None:
        raise RuntimeError("CpGPT human dependencies not resolved from cache")

    checkpoint_path = Path(resources.checkpoint_path)
    config_path = Path(resources.config_path)
    dependencies_path = Path(resources.dependencies_path)
    checkpoint_sha = sha256_file(checkpoint_path)
    configuration_sha = sha256_file(config_path)
    vocab_path = resources.vocab_path
    vocabulary_sha = sha256_file(Path(vocab_path)) if vocab_path is not None else None

    cpgpt_root = project_root / "vendor" / "cpgpt"
    source_commit = git_commit(cpgpt_root)

    loci = pd.read_parquet(
        resolved_loci,
        columns=["locus_id", "chromosome", "position", "canonical_key"],
    )
    location_index = _load_location_index(dependencies_path)
    loci_frame, mapped_source_rows = _build_loci_frame(loci, location_index=location_index)

    dna_mmap = _open_dna_mmap(
        dependencies_path,
        max_row_index=int(mapped_source_rows.max()) if mapped_source_rows.size else 0,
    )
    if mapped_source_rows.size:
        dna_rows = np.asarray(dna_mmap[mapped_source_rows], dtype=np.float32)
    else:
        dna_rows = np.empty((0, INPUT_DIMENSION), dtype=np.float32)

    if device.startswith("cuda") and not torch.cuda.is_available():
        device = "cpu"

    adapter = load_small_sequence_adapter(checkpoint_path, device=device)
    embeddings = _encode_batches(
        adapter,
        dna_rows,
        device=device,
        batch_size=batch_size,
    )
    validate_embeddings_array(embeddings, output_dimension=OUTPUT_DIMENSION)
    coverage = validate_loci_frame(loci_frame, n_mapped=int(embeddings.shape[0]))
    emb_stats = embedding_summary_stats(embeddings)

    paths = static_feature_store_paths(output_dir)
    paths.root.mkdir(parents=True, exist_ok=True)
    write_embeddings_zarr(paths.embeddings_path, embeddings)
    write_loci_index(paths.loci_path, loci_frame)

    command = export_command or (
        f"uv run --extra cpgpt mbs features export-cpgpt --feature-set-id {feature_set_id}"
    )
    manifest: dict[str, Any] = {
        "artifact_version": ARTIFACT_VERSION,
        "feature_set_id": feature_set_id,
        "source_model": DNA_LLM,
        "source_repository": SOURCE_REPOSITORY,
        "source_commit": source_commit,
        "checkpoint_sha256": checkpoint_sha,
        "configuration_sha256": configuration_sha,
        "vocabulary_sha256": vocabulary_sha,
        "context_length": CONTEXT_LENGTH,
        "genome_build": "GRCh38",
        "input_dimension": INPUT_DIMENSION,
        "output_dimension": OUTPUT_DIMENSION,
        "storage_dtype": STORAGE_DTYPE,
        "normalization": None,
        "n_loci": coverage["n_mapped"],
        "locus_table_sha256": loci_sha,
        "embedding_path": str(paths.embeddings_path),
        "locus_index_path": str(paths.loci_path),
        "export_command": command,
        "created_at": utc_now_iso(),
        "notes": (
            "CpGPT2M (small) sequence-adapter via dna_encoder / encode_sequence equivalent "
            "(mbs.static_features.cpgpt_adapter.SequenceAdapterMLP; avoids torchtune import); "
            "MBS 1-based cytosine positions mapped to CpGPT 0-based Ensembl keys; "
            "precomputed NTv2 human dependencies only (no DNA-LM regeneration)."
        ),
    }
    write_artifact(paths.artifact_path, manifest)

    stats = {
        **coverage,
        "mapping_rate": (coverage["n_mapped"] / coverage["n_loci"] if coverage["n_loci"] else 0.0),
        "embedding_stats": emb_stats,
        "checkpoint_path": str(checkpoint_path),
        "config_path": str(config_path),
        "dependencies_path": str(dependencies_path),
        "device": device,
    }
    if report_dir is not None:
        _write_report(report_dir, manifest=manifest, stats=stats)

    return CpGPTExportResult(
        feature_set_id=feature_set_id,
        output_dir=paths.root,
        report_dir=report_dir,
        manifest=manifest,
        stats=stats,
    )
