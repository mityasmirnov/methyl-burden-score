"""Virtual multi-store index and GSM overlap concordance for Hub pack matrices."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from mbs.annotation.manifest import write_json
from mbs.matrix.hub_pack import SUPPORTED_PACK_FAMILIES
from mbs.matrix.store import open_betas_zarr, read_sample_index

DEFAULT_OVERLAP_TOLERANCE = 1e-4
DEFAULT_COMPARE_LOCI = 256


@dataclass(frozen=True, slots=True)
class OverlapCheckResult:
    n_shared_gsm: int
    n_pairs_checked: int
    n_concordant: int
    n_discordant: int
    max_abs_diff: float
    discordant_pairs: tuple[tuple[str, str, str, float], ...]
    report: dict[str, Any]


def discover_full_pack_matrices(data_root: Path) -> list[tuple[str, str, Path]]:
    """Return ``(family, matrix_id, root)`` for existing ``matrix-hub-*-full-v1`` stores."""
    matrices_root = data_root / "canonical" / "matrices"
    found: list[tuple[str, str, Path]] = []
    for family in SUPPORTED_PACK_FAMILIES:
        matrix_id = f"matrix-hub-{family}-full-v1"
        root = matrices_root / matrix_id
        if (root / "betas.zarr").exists() and (root / "sample_index.parquet").is_file():
            found.append((family, matrix_id, root))
    return found


def build_hub_pack_matrix_index(
    data_root: Path,
    *,
    output_path: Path | None = None,
) -> pd.DataFrame:
    """Build virtual multi-store index (one row per family×GSM membership)."""
    rows: list[dict[str, Any]] = []
    for family, matrix_id, root in discover_full_pack_matrices(data_root):
        sample_index = read_sample_index(root / "sample_index.parquet")
        pheno_path = root / "sample_phenotypes.parquet"
        platform_by_id: dict[str, str] = {}
        if pheno_path.is_file():
            pheno = pd.read_parquet(pheno_path, columns=["sample_id", "platform"])
            for rec in pheno.drop_duplicates("sample_id").to_dict(orient="records"):
                platform_by_id[str(rec["sample_id"])] = str(rec.get("platform") or "")
        betas_path = str((root / "betas.zarr").resolve())
        for _, row in sample_index.iterrows():
            sid = str(row["sample_id"])
            rows.append(
                {
                    "family": family,
                    "matrix_id": matrix_id,
                    "sample_id": sid,
                    "row_index": int(row["row_index"]),
                    "platform": platform_by_id.get(sid, ""),
                    "betas_path": betas_path,
                }
            )
    frame = pd.DataFrame(rows)
    out = (
        output_path
        if output_path is not None
        else data_root / "canonical" / "matrices" / "hub_pack_matrix_index.parquet"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(out, index=False)
    return frame


def check_overlapping_gsm_betas(
    data_root: Path,
    *,
    index: pd.DataFrame | None = None,
    n_loci: int = DEFAULT_COMPARE_LOCI,
    tolerance: float = DEFAULT_OVERLAP_TOLERANCE,
    report_path: Path | None = None,
) -> OverlapCheckResult:
    """Compare betas for GSMs present in ≥2 pack matrices.

    Does **not** merge packs. Discordant pairs (max abs diff > tolerance) are
    reported; callers must not silently take the first pack.
    """
    if index is None:
        index_path = data_root / "canonical" / "matrices" / "hub_pack_matrix_index.parquet"
        if not index_path.is_file():
            index = build_hub_pack_matrix_index(data_root)
        else:
            index = pd.read_parquet(index_path)

    if index.empty:
        result = OverlapCheckResult(
            n_shared_gsm=0,
            n_pairs_checked=0,
            n_concordant=0,
            n_discordant=0,
            max_abs_diff=0.0,
            discordant_pairs=(),
            report={"status": "empty_index"},
        )
        if report_path is not None:
            write_json(report_path, result.report)
        return result

    counts = index.groupby("sample_id")["family"].nunique()
    shared_ids = counts[counts >= 2].index.astype(str).tolist()
    discordant: list[tuple[str, str, str, float]] = []
    max_abs = 0.0
    n_pairs = 0
    n_concordant = 0

    # Cache open arrays by matrix_id
    arrays: dict[str, Any] = {}
    for _, row in index.drop_duplicates("matrix_id").iterrows():
        mid = str(row["matrix_id"])
        arrays[mid] = open_betas_zarr(Path(str(row["betas_path"])))

    shared_set = set(shared_ids)
    by_sid: dict[str, list[dict[str, Any]]] = {}
    for rec in index.to_dict(orient="records"):
        sid = str(rec["sample_id"])
        if sid in shared_set:
            by_sid.setdefault(sid, []).append(rec)

    for sid in shared_ids:
        records = sorted(by_sid[sid], key=lambda r: str(r["family"]))
        for i in range(len(records)):
            for j in range(i + 1, len(records)):
                a = records[i]
                b = records[j]
                arr_a = arrays[str(a["matrix_id"])]
                arr_b = arrays[str(b["matrix_id"])]
                n_compare = min(n_loci, int(arr_a.shape[1]), int(arr_b.shape[1]))
                va = np.asarray(arr_a[int(a["row_index"]), :n_compare], dtype=np.float32)
                vb = np.asarray(arr_b[int(b["row_index"]), :n_compare], dtype=np.float32)
                both = np.isfinite(va) & np.isfinite(vb)
                if not both.any():
                    continue
                diff = float(np.max(np.abs(va[both] - vb[both])))
                max_abs = max(max_abs, diff)
                n_pairs += 1
                if diff > tolerance:
                    discordant.append((sid, str(a["family"]), str(b["family"]), diff))
                else:
                    n_concordant += 1

    report = {
        "n_shared_gsm": len(shared_ids),
        "n_pairs_checked": n_pairs,
        "n_concordant": n_concordant,
        "n_discordant": len(discordant),
        "max_abs_diff": max_abs,
        "tolerance": tolerance,
        "n_loci_compared": n_loci,
        "discordant_pairs": [
            {
                "sample_id": s,
                "family_a": fa,
                "family_b": fb,
                "max_abs_diff": d,
            }
            for s, fa, fb, d in discordant[:50]
        ],
        "status": "discordant" if discordant else "concordant",
        "merge_allowed": len(discordant) == 0,
    }
    if report_path is not None:
        write_json(report_path, report)
    return OverlapCheckResult(
        n_shared_gsm=len(shared_ids),
        n_pairs_checked=n_pairs,
        n_concordant=n_concordant,
        n_discordant=len(discordant),
        max_abs_diff=max_abs,
        discordant_pairs=tuple(discordant),
        report=report,
    )
