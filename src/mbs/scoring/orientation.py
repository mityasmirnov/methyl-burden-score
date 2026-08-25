"""ADR 0008 score-orientation anchor (predictive MBS, not constraint/LOEUF)."""

from __future__ import annotations

from typing import Any

import numpy as np
import torch
from torch import nn

from mbs.models import SeedMaskedLinearHead

ANCHOR_RECIPE = "pearson_vs_signed_gene_mean_mvalue"


def polarity_from_correlation(mbs: np.ndarray, signed_m: np.ndarray) -> str:
    """Return ``hyper_aligned`` or ``flipped`` from Pearson r of gene-mean scores."""
    a = np.asarray(mbs, dtype=np.float64).reshape(-1)
    b = np.asarray(signed_m, dtype=np.float64).reshape(-1)
    if a.shape != b.shape:
        raise ValueError(f"mbs and signed_m length mismatch: {a.shape} vs {b.shape}")
    if a.size < 2 or float(np.std(a)) == 0 or float(np.std(b)) == 0:
        return "hyper_aligned"
    r = float(np.corrcoef(a, b)[0, 1])
    return "flipped" if r < 0 else "hyper_aligned"


def signed_gene_mean_m(
    cpg_m: np.ndarray,
    cpg_to_gene: np.ndarray,
    n_genes: int,
    *,
    observed_mask: np.ndarray | None = None,
) -> np.ndarray:
    """Mean observed M-value per gene index from ragged CpG edges."""
    m = np.asarray(cpg_m, dtype=np.float64).reshape(-1)
    genes = np.asarray(cpg_to_gene, dtype=np.int64).reshape(-1)
    if m.shape != genes.shape:
        raise ValueError("cpg_m and cpg_to_gene length mismatch")
    if observed_mask is None:
        keep = np.isfinite(m)
    else:
        keep = np.asarray(observed_mask, dtype=bool).reshape(-1) & np.isfinite(m)
    sums = np.zeros(int(n_genes), dtype=np.float64)
    counts = np.zeros(int(n_genes), dtype=np.float64)
    if keep.any():
        g = genes[keep]
        vals = m[keep]
        np.add.at(sums, g, vals)
        np.add.at(counts, g, 1.0)
    out = np.zeros(int(n_genes), dtype=np.float64)
    nz = counts > 0
    out[nz] = sums[nz] / counts[nz]
    return out.astype(np.float32, copy=False)


def accumulate_signed_gene_mean_m(
    *,
    n_genes: int,
    cpg_m_batches: list[np.ndarray],
    cpg_to_gene_batches: list[np.ndarray],
) -> np.ndarray:
    """Pool signed gene-mean M across many samples' edge lists."""
    sums = np.zeros(int(n_genes), dtype=np.float64)
    counts = np.zeros(int(n_genes), dtype=np.float64)
    for cpg_m, cpg_to_gene in zip(cpg_m_batches, cpg_to_gene_batches, strict=True):
        m = np.asarray(cpg_m, dtype=np.float64).reshape(-1)
        genes = np.asarray(cpg_to_gene, dtype=np.int64).reshape(-1)
        keep = np.isfinite(m)
        if not keep.any():
            continue
        g = genes[keep]
        vals = m[keep]
        np.add.at(sums, g, vals)
        np.add.at(counts, g, 1.0)
    out = np.zeros(int(n_genes), dtype=np.float64)
    nz = counts > 0
    out[nz] = sums[nz] / counts[nz]
    return out.astype(np.float32, copy=False)


def gene_mean_mbs(mbs: np.ndarray, present: np.ndarray | None = None) -> np.ndarray:
    """Per-gene mean MBS over samples (optionally only where present)."""
    scores = np.asarray(mbs, dtype=np.float64)
    if scores.ndim == 1:
        return scores.astype(np.float32, copy=False)
    if present is None:
        return scores.mean(axis=0).astype(np.float32, copy=False)
    mask = np.asarray(present, dtype=bool)
    out = np.zeros(scores.shape[1], dtype=np.float64)
    for g in range(scores.shape[1]):
        col = scores[:, g]
        m = mask[:, g]
        out[g] = float(col[m].mean()) if m.any() else float(col.mean())
    return out.astype(np.float32, copy=False)


def apply_orientation(
    mbs: np.ndarray,
    *,
    signed_m: np.ndarray,
    head_weights: np.ndarray | None = None,
) -> dict[str, Any]:
    """Flip MBS (and optional head weights) so higher MBS tracks hypermethylation."""
    gene_mbs = gene_mean_mbs(mbs) if np.asarray(mbs).ndim == 2 else np.asarray(mbs)
    polarity = polarity_from_correlation(gene_mbs, signed_m)
    scores = np.asarray(mbs, dtype=np.float64)
    weights = None if head_weights is None else np.asarray(head_weights, dtype=np.float64).copy()
    if polarity == "flipped":
        scores = 1.0 - scores
        if weights is not None:
            weights = -weights
    return {
        "mbs": scores.astype(np.float32, copy=False),
        "head_weights": weights,
        "score_polarity": polarity,
        "anchor_recipe": ANCHOR_RECIPE,
        "score_family": "predictive_mbs",
    }


def orient_run_scores(
    mbs: np.ndarray,
    *,
    signed_m: np.ndarray,
    present: np.ndarray | None = None,
    head_weights: np.ndarray | None = None,
) -> dict[str, Any]:
    """Orient using per-gene mean MBS (present-aware) vs signed gene-mean M."""
    gene_mbs = gene_mean_mbs(mbs, present)
    return apply_orientation(gene_mbs, signed_m=signed_m, head_weights=head_weights)


def flip_phenotype_head_weights_(head: nn.Module) -> None:
    """In-place negate linear gene weights after a polarity flip (ADR 0008)."""
    if isinstance(head, SeedMaskedLinearHead):
        head.gene_weight.data.mul_(-1.0)
        return
    age_head = getattr(head, "age_head", None)
    tissue_head = getattr(head, "tissue_head", None)
    if isinstance(age_head, nn.Linear) and isinstance(tissue_head, SeedMaskedLinearHead):
        age_head.weight.data.mul_(-1.0)
        tissue_head.gene_weight.data.mul_(-1.0)
        for attr in ("sex_head", "disease_head", "cancer_head"):
            sub = getattr(head, attr, None)
            if isinstance(sub, nn.Linear):
                sub.weight.data.mul_(-1.0)
        return
    if isinstance(head, nn.Linear):
        head.weight.data.mul_(-1.0)
        return
    gw = getattr(head, "gene_weight", None)
    if isinstance(gw, torch.nn.Parameter):
        gw.data.mul_(-1.0)


def score_manifest(
    *,
    score_polarity: str,
    fold_id: str | None = None,
    restart_id: str | None = None,
) -> dict[str, Any]:
    return {
        "artifact_version": "1",
        "score_family": "predictive_mbs",
        "score_polarity": score_polarity,
        "anchor_recipe": ANCHOR_RECIPE,
        "fold_id": fold_id,
        "restart_id": restart_id,
        "notes": "Predictive sample x gene MBS; not a methylation constraint or LOEUF analogue.",
    }
