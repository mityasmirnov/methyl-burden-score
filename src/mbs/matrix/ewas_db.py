"""Read EWAS Data Hub All Data (``EWAS_db``) per-sample beta text files."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass(frozen=True, slots=True)
class EwasDbSampleFile:
    """One Hub sample beta file."""

    path: Path
    sample_id: str
    source_sample_id: str


@dataclass(frozen=True, slots=True)
class SampleBetaTable:
    """Parsed probe→beta table for one sample."""

    sample_id: str
    source_sample_id: str
    path: Path
    probe_ids: np.ndarray  # object / str
    betas: np.ndarray  # float64 before cast


def list_ewas_db_sample_files(source_dir: Path) -> list[EwasDbSampleFile]:
    """Return sorted ``*.txt`` sample files under an ``EWAS_db/{STUDY}/`` directory."""
    source_dir = source_dir.resolve()
    if not source_dir.is_dir():
        raise FileNotFoundError(f"EWAS_db study directory not found: {source_dir}")
    files = sorted(p for p in source_dir.iterdir() if p.is_file() and p.suffix == ".txt")
    if not files:
        raise FileNotFoundError(f"no *.txt sample files under {source_dir}")
    out: list[EwasDbSampleFile] = []
    for path in files:
        sample_id = path.stem
        if not sample_id:
            raise ValueError(f"empty sample id for {path}")
        out.append(
            EwasDbSampleFile(
                path=path,
                sample_id=sample_id,
                source_sample_id=sample_id,
            )
        )
    return out


def read_ewas_db_sample(path: Path, *, sample_id: str | None = None) -> SampleBetaTable:
    """Parse a Hub ``probe_id\\tbeta`` text file (no header).

    Finite betas outside ``[0, 1]`` are retained (never clipped) and counted by
    callers. Missing tokens are stored as NaN.
    """
    path = path.resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    sid = sample_id or path.stem
    frame = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["probe_id", "beta"],
        dtype={"probe_id": "string"},
        engine="c",
        na_values=["", "NA", "NaN", "nan", "."],
        keep_default_na=True,
    )
    if frame.empty:
        raise ValueError(f"empty sample beta file: {path}")
    probe_null = frame["probe_id"].isna()
    if bool(probe_null.any()):
        raise ValueError(f"null probe_id in {path}")
    probe_dup = frame["probe_id"].duplicated()
    if bool(probe_dup.any()):
        dup = str(frame.loc[probe_dup, "probe_id"].iloc[0])
        raise ValueError(f"duplicate probe_id {dup!r} in {path}")
    beta_series = pd.to_numeric(frame["beta"], errors="coerce")
    betas = np.asarray(beta_series, dtype=np.float64)
    probe_ids = frame["probe_id"].astype(str).to_numpy(dtype=object, copy=False)
    return SampleBetaTable(
        sample_id=sid,
        source_sample_id=sid,
        path=path,
        probe_ids=probe_ids,
        betas=betas,
    )


def beta_qc_stats(betas: np.ndarray) -> dict[str, float | int]:
    """Summarize beta values without modifying them."""
    finite = np.isfinite(betas)
    n_total = int(betas.size)
    n_missing = int((~finite).sum())
    if not finite.any():
        return {
            "n_values": n_total,
            "n_missing": n_missing,
            "n_finite": 0,
            "n_out_of_range": 0,
            "min": float("nan"),
            "max": float("nan"),
            "mean": float("nan"),
        }
    values = betas[finite]
    oor = (values < 0.0) | (values > 1.0)
    return {
        "n_values": n_total,
        "n_missing": n_missing,
        "n_finite": int(finite.sum()),
        "n_out_of_range": int(oor.sum()),
        "min": float(values.min()),
        "max": float(values.max()),
        "mean": float(values.mean()),
    }
