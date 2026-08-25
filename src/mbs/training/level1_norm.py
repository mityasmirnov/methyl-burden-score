"""Milestone 7D Level-1 fold-fitted MAD robust-z on train-fold M-values."""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Sequence
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from mbs.annotation.manifest import write_json

# Consistency factor: MAD -> approx. Gaussian sigma (STRATEGIC_PLAN / post-v0).
MAD_SCALE = 1.4826
_DEFAULT_SIGMA_MIN = 1e-6
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_FORMULA = "z=(M-median)/max(1.4826*MAD,sigma_min); train-fold only"


@dataclass(frozen=True, slots=True)
class Level1NormParams:
    """Per-locus robust reference fitted on the training fold."""

    mu: NDArray[np.float64]
    sigma: NDArray[np.float64]
    estimated: NDArray[np.bool_]
    locus_ids: NDArray[np.int64]
    sigma_min: float
    n_train_samples: int
    epsilon: float = 0.001
    fold_id: str | None = None
    run_id: str | None = None

    @property
    def n_loci(self) -> int:
        return int(self.mu.shape[0])

    @property
    def n_estimated(self) -> int:
        return int(self.estimated.sum())

    @property
    def n_unestimated(self) -> int:
        return int(self.n_loci - self.n_estimated)


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _sha256_array(arr: np.ndarray) -> str:
    return _sha256_bytes(np.ascontiguousarray(arr).tobytes())


def _as_2d(m: ArrayLike) -> NDArray[np.float64]:
    arr = np.asarray(m, dtype=np.float64)
    if arr.ndim != 2:
        raise ValueError(f"M matrix must be 2-D [n_samples, n_loci], got shape {arr.shape}")
    return arr


def fit_level1_from_betas(
    betas: ArrayLike,
    train_rows: Sequence[int],
    *,
    epsilon: float = 0.001,
    sigma_min: float = _DEFAULT_SIGMA_MIN,
    chunk_cols: int = 4096,
    locus_ids: ArrayLike | None = None,
    n_loci: int | None = None,
    fold_id: str | None = None,
    run_id: str | None = None,
) -> Level1NormParams:
    """Fit Level-1 from a beta matrix using only ``train_rows`` (column-chunked).

    Never materializes the full sample x locus dense block; column chunks are
    converted to M then summarized (ponytail: ``chunk_cols`` is the RAM ceiling).
    """
    # Deferred to avoid features <-> level1_norm import cycle.
    from mbs.training.features import beta_to_m_value  # noqa: PLC0415

    if sigma_min <= 0:
        raise ValueError("sigma_min must be positive")
    if chunk_cols < 1:
        raise ValueError("chunk_cols must be >= 1")
    rows = np.asarray(list(train_rows), dtype=np.int64)
    if rows.size == 0:
        raise ValueError("cannot fit Level-1 without train rows")
    full_loci = int(betas.shape[1])  # type: ignore[union-attr]
    if n_loci is None:
        n_loci = full_loci
    elif n_loci < 1:
        raise ValueError("n_loci must be >= 1")
    elif n_loci > full_loci:
        raise ValueError(f"n_loci={n_loci} exceeds beta matrix width {full_loci}")
    n_samples = int(rows.shape[0])

    if locus_ids is None:
        ids = np.arange(n_loci, dtype=np.int64)
    else:
        ids = np.asarray(locus_ids, dtype=np.int64).reshape(-1)
        if ids.shape[0] != n_loci:
            raise ValueError("locus_ids length must equal n_loci")

    mu = np.full(n_loci, np.nan, dtype=np.float64)
    sigma = np.full(n_loci, np.nan, dtype=np.float64)
    estimated = np.zeros(n_loci, dtype=bool)

    for start in range(0, n_loci, chunk_cols):
        stop = min(start + chunk_cols, n_loci)
        # Row-subset x column-chunk only.
        try:
            block_b = np.asarray(betas[np.ix_(rows, np.arange(start, stop))], dtype=np.float64)
        except Exception:
            # Fallback for stores that lack fancy indexing.
            block_b = np.asarray(
                [np.asarray(betas[int(r), start:stop], dtype=np.float64) for r in rows],
                dtype=np.float64,
            )
        block_m = beta_to_m_value(block_b, epsilon=epsilon).astype(np.float64, copy=False)
        obs = np.isfinite(block_m)
        width = stop - start
        for j in range(width):
            vals = block_m[obs[:, j], j]
            if vals.size == 0:
                continue
            med = float(np.median(vals))
            mad = float(np.median(np.abs(vals - med)))
            mu[start + j] = med
            sigma[start + j] = max(MAD_SCALE * mad, float(sigma_min))
            estimated[start + j] = True

    return Level1NormParams(
        mu=mu,
        sigma=sigma,
        estimated=estimated,
        locus_ids=ids,
        sigma_min=float(sigma_min),
        n_train_samples=n_samples,
        epsilon=float(epsilon),
        fold_id=fold_id,
        run_id=run_id,
    )


def fit_level1_robust_z(
    m: ArrayLike,
    *,
    observed: ArrayLike | None = None,
    sigma_min: float = _DEFAULT_SIGMA_MIN,
    chunk_cols: int = 4096,
    locus_ids: ArrayLike | None = None,
    epsilon: float = 0.001,
    fold_id: str | None = None,
    run_id: str | None = None,
) -> Level1NormParams:
    """Fit per-locus median and 1.4826×MAD on train-fold M-values.

    ``m`` is array-like ``[n_samples, n_loci]`` (numpy or Zarr-indexable).
    Columns are processed in chunks so Hub-scale fits never dense-stack the full
    matrix into RAM at once (ponytail: chunk width is the ceiling; raise
    ``chunk_cols`` if profiling shows I/O-bound).
    """
    if sigma_min <= 0:
        raise ValueError("sigma_min must be positive")
    if chunk_cols < 1:
        raise ValueError("chunk_cols must be >= 1")
    m_arr = _as_2d(m)
    n_samples, n_loci = m_arr.shape
    if observed is None:
        obs = np.isfinite(m_arr)
    else:
        obs = np.asarray(observed, dtype=bool)
        if obs.shape != m_arr.shape:
            raise ValueError("observed shape must match M")
        obs = obs & np.isfinite(m_arr)

    if locus_ids is None:
        ids = np.arange(n_loci, dtype=np.int64)
    else:
        ids = np.asarray(locus_ids, dtype=np.int64).reshape(-1)
        if ids.shape[0] != n_loci:
            raise ValueError("locus_ids length must equal n_loci")

    mu = np.full(n_loci, np.nan, dtype=np.float64)
    sigma = np.full(n_loci, np.nan, dtype=np.float64)
    estimated = np.zeros(n_loci, dtype=bool)

    for start in range(0, n_loci, chunk_cols):
        stop = min(start + chunk_cols, n_loci)
        # Materialize one column chunk only.
        block = np.asarray(m_arr[:, start:stop], dtype=np.float64)
        obs_block = np.asarray(obs[:, start:stop], dtype=bool)
        width = stop - start
        for j in range(width):
            vals = block[obs_block[:, j], j]
            if vals.size == 0:
                continue
            med = float(np.median(vals))
            mad = float(np.median(np.abs(vals - med)))
            mu[start + j] = med
            sigma[start + j] = max(MAD_SCALE * mad, float(sigma_min))
            estimated[start + j] = True

    return Level1NormParams(
        mu=mu,
        sigma=sigma,
        estimated=estimated,
        locus_ids=ids,
        sigma_min=float(sigma_min),
        n_train_samples=int(n_samples),
        epsilon=float(epsilon),
        fold_id=fold_id,
        run_id=run_id,
    )


def apply_level1_robust_z(
    m_row: ArrayLike,
    params: Level1NormParams,
    *,
    col_index: ArrayLike | None = None,
) -> tuple[NDArray[np.float32], NDArray[np.bool_]]:
    """Apply fold params to one sample's M row (or selected columns).

    Unestimated / novel loci → ``z=0``, ``norm_present=False`` (keep the CpG).
    """
    m = np.asarray(m_row, dtype=np.float64).reshape(-1)
    if col_index is None:
        cols = np.arange(m.shape[0], dtype=np.int64)
        if m.shape[0] != params.n_loci:
            raise ValueError(
                f"m_row length {m.shape[0]} != n_loci {params.n_loci} "
                "(pass col_index for ragged edges)"
            )
    else:
        cols = np.asarray(col_index, dtype=np.int64).reshape(-1)
        if cols.shape[0] != m.shape[0]:
            raise ValueError("col_index length must match m_row")
        if cols.size and (int(cols.min()) < 0 or int(cols.max()) >= params.n_loci):
            raise ValueError("col_index out of range for Level-1 params")

    z = np.zeros(m.shape[0], dtype=np.float32)
    present = np.zeros(m.shape[0], dtype=bool)
    finite = np.isfinite(m)
    for i in range(m.shape[0]):
        if not finite[i]:
            continue
        c = int(cols[i])
        if not params.estimated[c]:
            continue
        z[i] = np.float32((m[i] - params.mu[c]) / params.sigma[c])
        present[i] = True
    return z, present


def validate_fold_norm_manifest(manifest: dict[str, Any]) -> None:
    required = [
        "artifact_version",
        "formula",
        "epsilon",
        "sigma_min",
        "n_train_samples",
        "n_loci",
        "n_estimated",
        "n_unestimated",
        "mu_sha256",
        "sigma_sha256",
        "locus_ids_sha256",
        "mad_scale",
    ]
    missing = [k for k in required if k not in manifest]
    if missing:
        raise ValueError(f"fold_norm manifest missing keys: {missing}")
    for key in ("mu_sha256", "sigma_sha256", "locus_ids_sha256"):
        if not _SHA256_RE.fullmatch(str(manifest[key])):
            raise ValueError(f"{key} must be 64-char lowercase hex")


def persist_level1(
    run_root: Path,
    params: Level1NormParams,
    *,
    artifact_version: str = "fold-norm-level1-v1",
) -> dict[str, Any]:
    """Write ``fold_norm/`` under a run directory; never under canonical matrices."""
    out = Path(run_root) / "fold_norm"
    out.mkdir(parents=True, exist_ok=True)
    mu_path = out / "mu.npy"
    sigma_path = out / "sigma.npy"
    locus_path = out / "locus_ids.npy"
    np.save(mu_path, params.mu)
    np.save(sigma_path, params.sigma)
    np.save(locus_path, params.locus_ids)
    manifest: dict[str, Any] = {
        "artifact_version": artifact_version,
        "formula": _FORMULA,
        "fold_id": params.fold_id,
        "run_id": params.run_id,
        "epsilon": params.epsilon,
        "sigma_min": params.sigma_min,
        "mad_scale": MAD_SCALE,
        "n_train_samples": params.n_train_samples,
        "n_loci": params.n_loci,
        "n_estimated": params.n_estimated,
        "n_unestimated": params.n_unestimated,
        "mu_sha256": _sha256_array(params.mu),
        "sigma_sha256": _sha256_array(params.sigma),
        "locus_ids_sha256": _sha256_array(params.locus_ids),
        "mu_path": "fold_norm/mu.npy",
        "sigma_path": "fold_norm/sigma.npy",
        "locus_ids_path": "fold_norm/locus_ids.npy",
        "created_at": datetime.now(UTC).isoformat(),
    }
    validate_fold_norm_manifest(manifest)
    write_json(out / "manifest.json", manifest)
    return manifest


def load_level1(run_root: Path) -> tuple[Level1NormParams, dict[str, Any]]:
    """Load fold_norm artifacts and verify content hashes."""
    out = Path(run_root) / "fold_norm"
    manifest_path = out / "manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"missing fold_norm manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    validate_fold_norm_manifest(manifest)
    mu = np.load(out / "mu.npy")
    sigma = np.load(out / "sigma.npy")
    locus_ids = np.load(out / "locus_ids.npy")
    if _sha256_array(mu) != manifest["mu_sha256"]:
        raise ValueError("mu.npy sha256 mismatch vs fold_norm manifest")
    if _sha256_array(sigma) != manifest["sigma_sha256"]:
        raise ValueError("sigma.npy sha256 mismatch vs fold_norm manifest")
    if _sha256_array(locus_ids) != manifest["locus_ids_sha256"]:
        raise ValueError("locus_ids.npy sha256 mismatch vs fold_norm manifest")
    estimated = np.isfinite(mu) & np.isfinite(sigma)
    params = Level1NormParams(
        mu=np.asarray(mu, dtype=np.float64),
        sigma=np.asarray(sigma, dtype=np.float64),
        estimated=estimated,
        locus_ids=np.asarray(locus_ids, dtype=np.int64),
        sigma_min=float(manifest["sigma_min"]),
        n_train_samples=int(manifest["n_train_samples"]),
        epsilon=float(manifest["epsilon"]),
        fold_id=manifest.get("fold_id"),
        run_id=manifest.get("run_id"),
    )
    return params, manifest


def resolve_level1_config(config: dict[str, Any]) -> dict[str, Any]:
    """Read methylation Level-1 flags; fail loud on Level-2/3 stubs and mismatches."""
    feat = config.get("features", {}) if isinstance(config.get("features"), dict) else {}
    methyl = feat.get("methylation", {}) if isinstance(feat.get("methylation"), dict) else {}
    if "level2_probe_adapter" in methyl or methyl.get("level2_probe_adapter"):
        raise NotImplementedError(
            "Level-2 ProbeNormalizer is a documented ablation only (Milestone 7D); "
            "remove features.methylation.level2_probe_adapter"
        )
    if "level3_masked_ae" in methyl or methyl.get("level3_masked_ae"):
        raise NotImplementedError(
            "Level-3 masked AE is a documented ablation only (Milestone 7D); "
            "remove features.methylation.level3_masked_ae"
        )
    include_m = bool(methyl.get("m_value", True))
    robust = bool(methyl.get("robust_deviation", False))
    if robust and not include_m:
        raise ValueError(
            "features.methylation.robust_deviation=true requires m_value=true "
            "(Level-1 z is defined on M-values)"
        )
    return {
        "include_m_value": include_m,
        "include_robust_z": robust,
        "epsilon": float(methyl.get("epsilon", 0.001)),
        "sigma_min": float(methyl.get("sigma_min", _DEFAULT_SIGMA_MIN)),
    }
