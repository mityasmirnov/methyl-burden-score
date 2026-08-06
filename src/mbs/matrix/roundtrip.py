"""Round-trip verification from Hub raw text into the canonical matrix store."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from mbs.matrix.ewas_db import read_ewas_db_sample
from mbs.matrix.store import (
    MatrixStorePaths,
    matrix_store_paths,
    open_betas_zarr,
    read_locus_index,
    read_sample_index,
)


@dataclass(frozen=True, slots=True)
class RoundTripResult:
    ok: bool
    n_compared: int
    n_mismatch: int
    max_abs_diff: float
    sample_ids: tuple[str, ...]
    details: dict[str, Any]


def verify_roundtrip(
    source_dir: Path,
    output_dir: Path,
    *,
    sample_ids: list[str] | None = None,
    max_probes: int | None = 512,
    rtol: float = 0.0,
    atol: float = 0.0,
) -> RoundTripResult:
    """Compare stored float32 betas to values re-read from raw Hub files.

    Comparison uses float32 casting of the raw values (store dtype). Missing
    values (NaN) must match on both sides.
    """
    paths = matrix_store_paths(output_dir)
    sample_index = read_sample_index(paths.sample_index_path)
    locus_index = read_locus_index(paths.locus_index_path)
    betas = open_betas_zarr(paths.betas_path)

    if sample_ids is None:
        chosen = sample_index["sample_id"].astype(str).tolist()
    else:
        chosen = list(sample_ids)
        missing = set(chosen) - set(sample_index["sample_id"].astype(str))
        if missing:
            raise KeyError(f"sample_ids not in matrix: {sorted(missing)}")

    probe_ids = locus_index["probe_id"].astype(str).to_numpy()
    if max_probes is not None and max_probes < len(probe_ids):
        probe_ids = probe_ids[:max_probes]
        col_indices = np.arange(max_probes, dtype=np.int64)
    else:
        col_indices = locus_index["col_index"].to_numpy(dtype=np.int64)

    row_lookup = {
        str(sample_index.loc[i, "sample_id"]): int(sample_index.loc[i, "row_index"])
        for i in sample_index.index
    }

    n_compared = 0
    n_mismatch = 0
    max_abs_diff = 0.0
    per_sample: list[dict[str, Any]] = []

    for sample_id in chosen:
        raw_path = source_dir / f"{sample_id}.txt"
        table = read_ewas_db_sample(raw_path, sample_id=sample_id)
        raw_map = {
            str(pid): float(val) if np.isfinite(val) else float("nan")
            for pid, val in zip(table.probe_ids, table.betas, strict=True)
        }
        row = row_lookup[sample_id]
        stored = np.asarray(betas[row, col_indices], dtype=np.float32)
        expected_list: list[np.floating[Any]] = []
        for pid in probe_ids:
            if pid in raw_map:
                expected_list.append(np.float32(raw_map[pid]))
            else:
                expected_list.append(np.float32("nan"))
        expected = np.asarray(expected_list, dtype=np.float32)
        # Probes in the locus index must exist in the raw sample (union vocabulary
        # may include probes missing in a given sample → NaN).
        both_nan = np.isnan(stored) & np.isnan(expected)
        comparable = ~both_nan
        if comparable.any():
            abs_diff = np.abs(stored[comparable] - expected[comparable])
            mismatches = abs_diff > (atol + rtol * np.abs(expected[comparable]))
            n_mis = int(mismatches.sum())
            local_max = float(abs_diff.max()) if abs_diff.size else 0.0
        else:
            n_mis = 0
            local_max = 0.0
        # Also treat NaN asymmetry as mismatch
        nan_asym = int((np.isnan(stored) != np.isnan(expected)).sum())
        n_mis += nan_asym
        n_vals = len(probe_ids)
        n_compared += n_vals
        n_mismatch += n_mis
        max_abs_diff = max(max_abs_diff, local_max)
        per_sample.append(
            {
                "sample_id": sample_id,
                "n_compared": n_vals,
                "n_mismatch": n_mis,
                "max_abs_diff": local_max,
            }
        )

    return RoundTripResult(
        ok=n_mismatch == 0,
        n_compared=n_compared,
        n_mismatch=n_mismatch,
        max_abs_diff=max_abs_diff,
        sample_ids=tuple(chosen),
        details={"per_sample": per_sample, "paths": _paths_dict(paths)},
    )


def _paths_dict(paths: MatrixStorePaths) -> dict[str, str]:
    return {
        "root": str(paths.root),
        "betas_path": str(paths.betas_path),
        "sample_index_path": str(paths.sample_index_path),
        "locus_index_path": str(paths.locus_index_path),
        "manifest_path": str(paths.manifest_path),
    }
