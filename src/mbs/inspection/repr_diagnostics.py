"""Representation diagnostics for Stage A FlatDeepSet / ablation MBS matrices."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

SATURATION_LO = 0.05
SATURATION_HI = 0.95
CONST_SD_EPS = 1e-4


def compute_mbs_repr_diagnostics(
    mbs: np.ndarray,
    present: np.ndarray | None = None,
    *,
    sample_mean_m: np.ndarray | None = None,
    head_weight_l2: float | None = None,
    best_epoch: int | None = None,
    best_val_tissue_f1: float | None = None,
) -> dict[str, Any]:
    """Compute gene-score variance / saturation / constancy (+ optional anchors).

    Parameters
    ----------
    mbs:
        ``(n_samples, n_genes)`` MBS scores (raw encoder export preferred).
    present:
        Optional bool/uint8 mask, same shape. Missing entries are ignored.
    sample_mean_m:
        Optional ``(n_samples,)`` mean M-value over the gene-linked CpG panel.
        When provided, ``corr_mean_m`` is Pearson r between per-sample mean MBS
        and this vector (finite pairs only).
    """
    scores = np.asarray(mbs, dtype=np.float64)
    if scores.ndim != 2:
        raise ValueError(f"mbs must be 2D, got shape {scores.shape}")
    if present is None:
        mask = np.isfinite(scores)
    else:
        mask = np.asarray(present, dtype=bool) & np.isfinite(scores)

    vals = scores[mask]
    gene_score_sd = float(np.std(vals)) if vals.size else float("nan")
    saturation_fraction = (
        float(((vals <= SATURATION_LO) | (vals >= SATURATION_HI)).mean()) if vals.size else float("nan")
    )

    masked = np.where(mask, scores, np.nan)
    gene_sds = np.nanstd(masked, axis=0)
    gene_sds = np.where(np.isfinite(gene_sds), gene_sds, 0.0)
    const_score_fraction = float((gene_sds < CONST_SD_EPS).mean()) if gene_sds.size else float("nan")
    mean_per_gene_sd = float(np.mean(gene_sds)) if gene_sds.size else float("nan")
    n_genes = int(scores.shape[1])

    corr_mean_m = float("nan")
    if sample_mean_m is not None:
        sm = np.asarray(sample_mean_m, dtype=np.float64).reshape(-1)
        if sm.shape[0] != scores.shape[0]:
            raise ValueError(
                f"sample_mean_m length {sm.shape[0]} != n_samples {scores.shape[0]}"
            )
        sample_mean_mbs = np.nanmean(masked, axis=1)
        ok = np.isfinite(sample_mean_mbs) & np.isfinite(sm)
        if int(ok.sum()) >= 2 and float(np.std(sample_mean_mbs[ok])) > 0 and float(np.std(sm[ok])) > 0:
            corr_mean_m = float(np.corrcoef(sample_mean_mbs[ok], sm[ok])[0, 1])

    out: dict[str, Any] = {
        "gene_score_sd": gene_score_sd,
        "mean_per_gene_sd": mean_per_gene_sd,
        "saturation_fraction": saturation_fraction,
        "constant_score_fraction": const_score_fraction,
        "corr_mean_m": corr_mean_m,
        "n_scores": int(vals.size),
        "n_genes": int(n_genes),
        "n_samples": int(scores.shape[0]),
        "score_min": float(vals.min()) if vals.size else float("nan"),
        "score_max": float(vals.max()) if vals.size else float("nan"),
        "score_mean": float(vals.mean()) if vals.size else float("nan"),
    }
    if head_weight_l2 is not None:
        out["head_weight_l2"] = float(head_weight_l2)
    if best_epoch is not None:
        out["best_epoch"] = int(best_epoch)
    if best_val_tissue_f1 is not None:
        out["best_val_tissue_f1"] = float(best_val_tissue_f1)
    return out


def head_weight_l2_from_checkpoint(ckpt_path: Path | str) -> float | None:
    """L2 norm of all tensors in ``head_state`` (or None if missing)."""
    path = Path(ckpt_path)
    if not path.is_file():
        return None
    import torch

    blob = torch.load(path, map_location="cpu", weights_only=False)
    if not isinstance(blob, dict):
        return None
    head = blob.get("head_state")
    if not isinstance(head, dict) or not head:
        return None
    total = 0.0
    for t in head.values():
        total += float(np.sum(np.asarray(t.detach().cpu().numpy(), dtype=np.float64) ** 2))
    return float(np.sqrt(total))


def cache_sample_mean_m_gene_panel(
    *,
    betas_zarr: Path | str,
    gene_col_indices: np.ndarray,
    out_path: Path | str,
    chunk_rows: int = 256,
) -> np.ndarray:
    """Mean M-value over gene-linked CpG columns per sample; cache to ``out_path``."""
    import zarr

    out = Path(out_path)
    if out.is_file():
        return np.load(out)

    cols = np.asarray(gene_col_indices, dtype=np.int64)
    z = zarr.open(str(betas_zarr), mode="r")
    n_samples = int(z.shape[0])
    means = np.zeros(n_samples, dtype=np.float64)
    for start in range(0, n_samples, chunk_rows):
        stop = min(start + chunk_rows, n_samples)
        block = np.asarray(z.oindex[start:stop, cols], dtype=np.float64)
        # beta -> M = logit(beta); clip for numerical stability
        b = np.clip(block, 1e-6, 1.0 - 1e-6)
        mvals = np.log(b / (1.0 - b))
        means[start:stop] = np.nanmean(mvals, axis=1)
    out.parent.mkdir(parents=True, exist_ok=True)
    np.save(out, means.astype(np.float32))
    return means.astype(np.float32, copy=False)
