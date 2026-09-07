"""Virtual Hub nine-pack multi-store (Milestone 7E′).

Routes each GSM to one pack Zarr via a fixed priority order. Canonical locus
columns are the intersection of pack ``locus_id``s, ordered as the frozen age
pack. No dense union Zarr is written.
"""

from __future__ import annotations

import json
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from mbs.annotation.manifest import git_commit, sha256_file, utc_now_iso, write_json
from mbs.matrix.hub_pack_index import build_hub_pack_matrix_index, discover_full_pack_matrices
from mbs.matrix.store import (
    ARTIFACT_VERSION,
    DEFAULT_DTYPE,
    MATRIX_ORIENTATION,
    MISSING_VALUE_ENCODING,
    matrix_store_paths,
    open_betas_zarr,
    read_locus_index,
    write_locus_index,
    write_sample_index,
)
from mbs.platform_id import normalize_platform

VIRTUAL_MATRIX_ID = "matrix-hub-nine-pack-virtual-v1"
VIRTUAL_KIND = "virtual_multi_store"
CANONICAL_AGE_MATRIX_ID = "matrix-hub-age-full-v1"

# Overlap betas are concordant; prefer packs that match the 5d ATS first-seen order.
PACK_PRIORITY: tuple[str, ...] = (
    "age",
    "tissue",
    "sex",
    "disease",
    "cancer",
    "blood",
    "brain",
    "bmi",
    "ancestry",
)
_PRIORITY_RANK = {fam: i for i, fam in enumerate(PACK_PRIORITY)}


@dataclass(frozen=True, slots=True)
class VirtualHubBuildResult:
    matrix_id: str
    output_dir: Path
    n_samples: int
    n_loci: int
    n_source_matrices: int
    route_path: Path
    stats: dict[str, Any]


class RoutedBetas:
    """Zarr-like ``[n_samples, n_loci]`` view over multiple pack stores.

    Supports ``arr[row]``, ``arr[row, :]``, ``arr[row, start:stop]``,
    ``arr[:, :n]`` / ``arr[row_slice, col_slice]`` (dense training idiom), and
    ``arr[np.ix_(rows, cols)]`` for Level-1 fitting.
    """

    def __init__(
        self,
        *,
        route: pd.DataFrame,
        pack_arrays: dict[str, Any],
        pack_col_maps: dict[str, np.ndarray],
        n_loci: int,
    ) -> None:
        if route.empty:
            raise ValueError("route table is empty")
        self._route = route.sort_values("row_index").reset_index(drop=True)
        self._pack_arrays = pack_arrays
        self._pack_col_maps = pack_col_maps
        self._n_samples = len(self._route)
        self._n_loci = int(n_loci)
        max_row = int(self._route["row_index"].to_numpy().max())
        if self._n_samples != max_row + 1:
            raise ValueError("route row_index must be dense 0..n-1")

    @property
    def shape(self) -> tuple[int, int]:
        return (self._n_samples, self._n_loci)

    def _normalize_row_ids(self, rows: Any) -> np.ndarray:
        if isinstance(rows, slice):
            start = 0 if rows.start is None else int(rows.start)
            stop = self._n_samples if rows.stop is None else int(rows.stop)
            step = 1 if rows.step is None else int(rows.step)
            return np.arange(start, stop, step, dtype=np.int64)
        if isinstance(rows, (int, np.integer)):
            return np.asarray([int(rows)], dtype=np.int64)
        if isinstance(rows, np.ndarray):
            if rows.ndim == 2:
                return np.asarray(rows[:, 0], dtype=np.int64)
            return np.asarray(rows, dtype=np.int64).reshape(-1)
        return np.asarray(rows, dtype=np.int64).reshape(-1)

    def _normalize_col_ids(self, cols: Any) -> np.ndarray:
        if isinstance(cols, slice):
            start = 0 if cols.start is None else int(cols.start)
            stop = self._n_loci if cols.stop is None else int(cols.stop)
            step = 1 if cols.step is None else int(cols.step)
            return np.arange(start, stop, step, dtype=np.int64)
        if isinstance(cols, (int, np.integer)):
            return np.asarray([int(cols)], dtype=np.int64)
        if isinstance(cols, np.ndarray):
            if cols.ndim == 2:
                return np.asarray(cols[0, :], dtype=np.int64)
            return np.asarray(cols, dtype=np.int64).reshape(-1)
        return np.asarray(cols, dtype=np.int64).reshape(-1)

    def _block(self, row_ids: np.ndarray, col_ids: np.ndarray) -> np.ndarray:
        """Materialize ``(n_rows, n_cols)`` via pack-batched Zarr reads.

        Rows are processed in chunks per pack so cascade ``[:, :max_loci]`` loads
        (tens of thousands of samples × tens of thousands of columns) stay
        tractable. Contiguous virtual prefixes prefer full pack-row reads when
        that is cheaper than a huge fancy column index.
        """
        row_ids = np.asarray(row_ids, dtype=np.int64).reshape(-1)
        col_ids = np.asarray(col_ids, dtype=np.int64).reshape(-1)
        block = np.empty((row_ids.shape[0], col_ids.shape[0]), dtype=np.float32)
        if row_ids.size == 0 or col_ids.size == 0:
            return block
        contiguous_prefix = bool(
            col_ids.size > 0
            and int(col_ids[0]) == 0
            and np.array_equal(col_ids, np.arange(col_ids.size, dtype=np.int64))
        )
        # Dense training prefixes ([:, :max_loci]) and near-full-width requests:
        # read full pack rows then gather — fewer Zarr seeks than fancy 65k cols.
        wide = contiguous_prefix or col_ids.shape[0] >= max(self._n_loci // 2, 1)
        by_pack: dict[str, list[tuple[int, int]]] = {}
        for i, r in enumerate(row_ids.tolist()):
            rec = self._route.iloc[int(r)]
            mid = str(rec["matrix_id"])
            src_row = int(rec["src_row_index"])
            by_pack.setdefault(mid, []).append((i, src_row))
        # Bound peak RAM: full pack row is ~n_pack_loci floats.
        row_chunk = 128 if wide else (64 if col_ids.shape[0] >= 8192 else 256)
        n_packs = len(by_pack)
        for pack_i, (mid, pairs) in enumerate(by_pack.items()):
            if row_ids.shape[0] >= 4096 and col_ids.shape[0] >= 8192:
                print(
                    f"[RoutedBetas] pack {pack_i + 1}/{n_packs} {mid} "
                    f"rows={len(pairs)} cols={col_ids.shape[0]} "
                    f"chunk={row_chunk} wide={wide}",
                    flush=True,
                )
            arr = self._pack_arrays[mid]
            col_map = self._pack_col_maps[mid]
            src_cols = col_map[col_ids]
            valid = src_cols >= 0
            out_idx = np.where(valid)[0]
            valid_src = src_cols[valid]
            if not valid.any():
                order = [i for i, _ in pairs]
                block[order, :] = np.nan
                continue
            for start in range(0, len(pairs), row_chunk):
                chunk = pairs[start : start + row_chunk]
                order = [i for i, _ in chunk]
                src_rows = [sr for _, sr in chunk]
                if wide:
                    dense = np.asarray(arr[src_rows], dtype=np.float32)
                    gathered = dense[:, valid_src]
                else:
                    gathered = np.asarray(
                        arr[np.ix_(src_rows, valid_src.tolist())], dtype=np.float32
                    )
                out = np.full((len(order), col_ids.shape[0]), np.nan, dtype=np.float32)
                out[:, out_idx] = gathered
                block[order, :] = out
        return block

    def _row_vector(self, row: int, col_sl: slice | None = None) -> np.ndarray:
        col_ids = (
            np.arange(self._n_loci, dtype=np.int64)
            if col_sl is None
            else self._normalize_col_ids(col_sl)
        )
        return self._block(np.asarray([int(row)], dtype=np.int64), col_ids)[0]

    def __getitem__(self, key: Any) -> Any:
        if isinstance(key, tuple) and len(key) == 2:
            rows, cols = key
            # Dense training idiom: arr[:, :n] / arr[start:stop, :n] / np.ix_.
            if isinstance(rows, slice) or isinstance(rows, np.ndarray):
                row_ids = self._normalize_row_ids(rows)
                col_ids = self._normalize_col_ids(cols)
                return self._block(row_ids, col_ids)
            if isinstance(rows, (int, np.integer)):
                if isinstance(cols, slice):
                    return self._row_vector(int(rows), cols)
                if isinstance(cols, (int, np.integer)):
                    return float(self._row_vector(int(rows))[int(cols)])
                col_ids = self._normalize_col_ids(cols)
                return self._block(np.asarray([int(rows)], dtype=np.int64), col_ids)[0]
            raise TypeError(f"unsupported RoutedBetas index: {type(rows)}, {type(cols)}")
        if isinstance(key, (int, np.integer)):
            return self._row_vector(int(key))
        if isinstance(key, slice):
            row_ids = self._normalize_row_ids(key)
            col_ids = np.arange(self._n_loci, dtype=np.int64)
            return self._block(row_ids, col_ids)
        raise TypeError(f"unsupported RoutedBetas key: {type(key)}")


def is_virtual_matrix_manifest(manifest: dict[str, Any]) -> bool:
    return str(manifest.get("kind", "")) == VIRTUAL_KIND


def read_matrix_manifest(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def open_betas_for_matrix(root: Path) -> Any:
    """Open a dense Zarr or a virtual RoutedBetas store."""
    paths = matrix_store_paths(root)
    if not paths.manifest_path.is_file():
        return open_betas_zarr(paths.betas_path)
    manifest = read_matrix_manifest(paths.manifest_path)
    if not is_virtual_matrix_manifest(manifest):
        return open_betas_zarr(paths.betas_path)
    return open_routed_betas(root, manifest=manifest)


def open_routed_betas(root: Path, *, manifest: dict[str, Any] | None = None) -> RoutedBetas:
    paths = matrix_store_paths(root)
    man = manifest if manifest is not None else read_matrix_manifest(paths.manifest_path)
    route_rel = man.get("route_path")
    route_path = Path(str(route_rel)) if route_rel else root / "route.parquet"
    if not route_path.is_file():
        raise FileNotFoundError(f"virtual route missing: {route_path}")
    route = pd.read_parquet(route_path)
    locus = read_locus_index(paths.locus_index_path)
    canonical_ids = locus["locus_id"].to_numpy()
    n_loci = len(canonical_ids)

    pack_arrays: dict[str, Any] = {}
    pack_col_maps: dict[str, np.ndarray] = {}
    for mid in sorted({str(x) for x in route["matrix_id"].tolist()}):
        # Prefer path from route; fall back to canonical layout.
        sub = route.loc[route["matrix_id"].astype(str) == mid].iloc[0]
        betas_path = Path(str(sub["betas_path"]))
        pack_root = betas_path.parent if betas_path.name == "betas.zarr" else betas_path
        if (pack_root / "betas.zarr").exists():
            pack_arrays[mid] = open_betas_zarr(pack_root / "betas.zarr")
            pack_locus = read_locus_index(pack_root / "locus_index.parquet")
        else:
            pack_arrays[mid] = open_betas_zarr(betas_path)
            pack_locus = read_locus_index(betas_path.parent / "locus_index.parquet")
        id_to_col = {
            int(lid): int(col)
            for lid, col in zip(
                pack_locus["locus_id"].astype(np.int64),
                pack_locus["col_index"].astype(np.int64),
                strict=True,
            )
        }
        col_map = np.full(n_loci, -1, dtype=np.int64)
        for i, lid in enumerate(canonical_ids.astype(np.int64)):
            col_map[i] = id_to_col.get(int(lid), -1)
        pack_col_maps[mid] = col_map

    return RoutedBetas(
        route=route,
        pack_arrays=pack_arrays,
        pack_col_maps=pack_col_maps,
        n_loci=n_loci,
    )


def _choose_canonical_locus(
    data_root: Path,
    discovered: list[tuple[str, str, Path]],
) -> pd.DataFrame:
    by_id = {mid: root for _fam, mid, root in discovered}
    age_root = by_id.get(CANONICAL_AGE_MATRIX_ID)
    if age_root is None:
        # Fall back to highest-priority discovered pack.
        for fam in PACK_PRIORITY:
            for f, _mid, root in discovered:
                if f == fam:
                    age_root = root
                    break
            if age_root is not None:
                break
    if age_root is None:
        raise FileNotFoundError("no Hub full-pack matrices found for locus order")
    age_locus = read_locus_index(age_root / "locus_index.parquet")
    # Intersection with every pack's locus set.
    keep = set(age_locus["locus_id"].astype(np.int64).tolist())
    for _fam, _mid, root in discovered:
        other = read_locus_index(root / "locus_index.parquet")
        keep &= set(other["locus_id"].astype(np.int64).tolist())
    if not keep:
        raise ValueError("empty locus intersection across Hub packs")
    mask = age_locus["locus_id"].astype(np.int64).isin(list(keep))
    return age_locus.loc[mask].reset_index(drop=True)


def build_virtual_hub_store(
    *,
    data_root: Path,
    output_matrix_id: str = VIRTUAL_MATRIX_ID,
    output_dir: Path | None = None,
    index: pd.DataFrame | None = None,
) -> VirtualHubBuildResult:
    """Build routing + indices for the Hub nine-pack virtual cohort."""
    data_root = data_root.resolve()
    discovered = discover_full_pack_matrices(data_root)
    if not discovered:
        raise FileNotFoundError(f"no matrix-hub-*-full-v1 under {data_root}")

    if index is None:
        index_path = data_root / "canonical" / "matrices" / "hub_pack_matrix_index.parquet"
        if index_path.is_file():
            index = pd.read_parquet(index_path)
        else:
            index = build_hub_pack_matrix_index(data_root)

    if index.empty:
        raise ValueError("hub pack matrix index is empty")

    work = index.copy()
    work["family"] = work["family"].astype(str)
    work["sample_id"] = work["sample_id"].astype(str)
    work["matrix_id"] = work["matrix_id"].astype(str)
    work["priority"] = work["family"].map(lambda f: _PRIORITY_RANK.get(str(f), 999))  # type: ignore[misc]
    work = work.sort_values(["sample_id", "priority", "family"]).drop_duplicates(
        "sample_id", keep="first"
    )
    work = work.sort_values(["sample_id"]).reset_index(drop=True)
    ordered_ids = work["sample_id"].astype(str).tolist()
    n_samples = len(ordered_ids)

    locus = _choose_canonical_locus(data_root, discovered)
    n_loci = len(locus)

    out_dir = (
        output_dir
        if output_dir is not None
        else data_root / "canonical" / "matrices" / output_matrix_id
    ).resolve()
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = matrix_store_paths(out_dir)

    write_sample_index(
        paths.sample_index_path,
        sample_ids=ordered_ids,
        source_sample_ids=ordered_ids,
    )
    status = (
        locus["annotation_status"].to_numpy(dtype=np.int8)
        if "annotation_status" in locus.columns
        else None
    )
    write_locus_index(
        paths.locus_index_path,
        locus_ids=locus["locus_id"].to_numpy(),
        canonical_keys=locus["canonical_key"].to_numpy(),
        probe_ids=locus["probe_id"].to_numpy(),
        annotation_status=status,
    )

    route_rows: list[dict[str, Any]] = []
    for row_i, rec in enumerate(work.to_dict(orient="records")):
        platform = normalize_platform(rec.get("platform")) or str(rec.get("platform") or "")
        route_rows.append(
            {
                "row_index": row_i,
                "sample_id": str(rec["sample_id"]),
                "family": str(rec["family"]),
                "matrix_id": str(rec["matrix_id"]),
                "src_row_index": int(rec["row_index"]),
                "betas_path": str(rec["betas_path"]),
                "platform": platform,
            }
        )
    route = pd.DataFrame(route_rows)
    route_path = out_dir / "route.parquet"
    route.to_parquet(route_path, index=False)

    source_mids = sorted({str(x) for x in route["matrix_id"].tolist()})
    source_files = []
    for mid in source_mids:
        src_root = data_root / "canonical" / "matrices" / mid
        manifest_p = src_root / "matrix_manifest.json"
        if manifest_p.is_file():
            source_files.append(
                {
                    "path": str(manifest_p.resolve()),
                    "sha256": sha256_file(manifest_p),
                    "byte_size": int(manifest_p.stat().st_size),
                }
            )
    # Always include the route itself.
    source_files.append(
        {
            "path": str(route_path.resolve()),
            "sha256": sha256_file(route_path),
            "byte_size": int(route_path.stat().st_size),
        }
    )

    family_counts = route["family"].value_counts().to_dict()
    manifest = {
        "artifact_version": ARTIFACT_VERSION,
        "matrix_id": output_matrix_id,
        "study_id": "HUB_NINE_PACK_VIRTUAL",
        "platform_id": "HM450",
        "processing_level": "gmqn_beta",
        "genome_build": "GRCh38",
        "shape": [n_samples, n_loci],
        "dtype": DEFAULT_DTYPE,
        "chunks": None,
        "compression": None,
        "missing_value_encoding": MISSING_VALUE_ENCODING,
        "matrix_orientation": MATRIX_ORIENTATION,
        "matrix_path": str(route_path.resolve()),
        "sample_index_path": str(paths.sample_index_path.resolve()),
        "locus_index_path": str(paths.locus_index_path.resolve()),
        "sample_index_sha256": sha256_file(paths.sample_index_path),
        "locus_index_sha256": sha256_file(paths.locus_index_path),
        "source_files": source_files,
        "conversion_commit": git_commit(Path(__file__).resolve().parents[2]),
        "created_at": utc_now_iso(),
        "kind": VIRTUAL_KIND,
        "route_path": str(route_path.resolve()),
        "source_matrix_ids": source_mids,
        "notes": (
            "Virtual multi-store: betas live in per-pack Zarrs; route.parquet "
            f"maps row→pack. Pack priority: {','.join(PACK_PRIORITY)}. "
            f"Family routing counts: {family_counts}."
        ),
    }
    # validate via store helper but allow kind/route extras — write JSON directly
    # after a soft required-field check (schema allows additional notes only;
    # kind/route_path stored and read by open_routed_betas).
    _validate_virtual_manifest(manifest)
    write_json(paths.manifest_path, manifest)

    stats = {
        "n_samples": n_samples,
        "n_loci": n_loci,
        "n_source_matrices": len(source_mids),
        "family_routing_counts": {str(k): int(v) for k, v in family_counts.items()},
        "pack_priority": list(PACK_PRIORITY),
    }
    return VirtualHubBuildResult(
        matrix_id=output_matrix_id,
        output_dir=out_dir,
        n_samples=n_samples,
        n_loci=n_loci,
        n_source_matrices=len(source_mids),
        route_path=route_path,
        stats=stats,
    )


def _validate_virtual_manifest(manifest: dict[str, Any]) -> None:
    """Required fields for virtual manifests (extends dense matrix rules lightly)."""
    required = [
        "artifact_version",
        "matrix_id",
        "study_id",
        "platform_id",
        "processing_level",
        "genome_build",
        "shape",
        "dtype",
        "missing_value_encoding",
        "matrix_path",
        "sample_index_path",
        "locus_index_path",
        "source_files",
        "conversion_commit",
        "created_at",
        "kind",
        "route_path",
    ]
    missing = [k for k in required if k not in manifest]
    if missing:
        raise ValueError(f"virtual matrix manifest missing keys: {missing}")
    if manifest["genome_build"] != "GRCh38":
        raise ValueError("genome_build must be GRCh38")
    if manifest["kind"] != VIRTUAL_KIND:
        raise ValueError(f"kind must be {VIRTUAL_KIND}")
    for key in ("matrix_path", "sample_index_path", "locus_index_path", "route_path"):
        if not str(manifest[key]).startswith("/data/"):
            raise ValueError(f"{key} must be an absolute /data path")
