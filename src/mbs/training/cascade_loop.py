"""Milestone 7F cascade trainer: RBS→gene MBS + orphan RBS + leftover direct."""

from __future__ import annotations

import hashlib
import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal, cast

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F

from mbs.annotation.manifest import write_json
from mbs.matrix.store import (
    matrix_store_paths,
    open_betas_zarr,
    read_locus_index,
    read_sample_index,
)
from mbs.models import CascadeDeepSet
from mbs.segment_ops import PoolName
from mbs.training.cascade_assign import (
    CascadeAssignment,
    assignment_gene_linked_only,
    build_cascade_assignment,
    gene_linked_col_index,
)
from mbs.training.cascade_scores import (
    FusionBlockMode,
    fusion_feature_matrix,
    load_cascade_score_blocks,
    write_cascade_score_dir,
    _write_array,
)
from mbs.training.checkpoint_selection import (
    validation_rank as _validation_rank,
    validation_rank_age_primary as _validation_rank_age_primary,
)
from mbs.training.dev_cv import DEFAULT_SPLIT_ID, load_frozen_folds
from mbs.training.direct_cpg import direct_cpg_design_matrix, fit_direct_elasticnet
from mbs.training.features import beta_to_m_value
from mbs.training.late_fusion import evaluate_late_fusion
from mbs.training.locus_gene import load_graph_tables
from mbs.training.loop import load_experiment_config, resolve_device
from mbs.training.multitask import MultitaskHeads
from mbs.training.phenotypes import load_multitask_phenotypes
from mbs.training.transparent_baselines import (
    evaluate_multitask_predictions,
    run_elasticnet_multitask,
)
from mbs.training.run_artifacts import run_dir

PrimaryEvaluation = Literal["late_fusion", "mbs_e2e"]

AGE_PRIMARY_SELECTION = "validation_age_mae_then_tissue_f1_then_sex_auroc"
TISSUE_PRIMARY_SELECTION = "validation_tissue_macro_f1_then_age_mae"


def _seed_mask_tensor(mask: np.ndarray | torch.Tensor | None) -> torch.Tensor | None:
    """Coerce a numpy/tensor seed mask to a CPU float32 tensor (validation in heads)."""
    if mask is None:
        return None
    if isinstance(mask, torch.Tensor):
        return mask.detach().to(torch.float32).cpu()
    return torch.as_tensor(np.asarray(mask), dtype=torch.float32)


def _seed_mask_meta(mask: torch.Tensor | None) -> dict[str, Any] | None:
    """Hash + per-row gene counts for a seed mask (auditability, ADR 0011)."""
    if mask is None:
        return None
    arr = np.ascontiguousarray(mask.detach().cpu().numpy().astype(np.float32))
    rows = arr if arr.ndim > 1 else arr.reshape(1, -1)
    return {
        "sha256": hashlib.sha256(arr.tobytes()).hexdigest(),
        "shape": list(arr.shape),
        "n_genes_selected": [int(x) for x in rows.sum(axis=1).tolist()],
    }


def _encoder_grad_norm(model: torch.nn.Module) -> float:
    """L2 norm of encoder parameter gradients (0 if none)."""
    total = 0.0
    for p in model.parameters():
        if p.grad is None:
            continue
        total += float(p.grad.detach().pow(2).sum().item())
    return float(total**0.5)


@dataclass(frozen=True, slots=True)
class CascadeTrainResult:
    metrics: dict[str, Any]
    score_dir: Path
    report_dir: Path


def _set_seed(seed: int) -> None:
    torch.manual_seed(seed)
    np.random.seed(seed)


def make_synthetic_cascade_tables(
    *,
    seed: int = 0,
) -> dict[str, Any]:
    """Tiny graph + dense betas for fixture tests (gene RBS, orphan RBS, leftover, fake TBS)."""
    rng = np.random.default_rng(seed)
    genes = pd.DataFrame(
        {
            "gene_id": ["ENSG1", "ENSG2"],
            "chromosome": ["chr1", "chr1"],
            "start": [100, 5000],
            "end": [400, 5500],
            "strand": ["+", "+"],
            "gene_name": ["G1", "G2"],
            "gene_type": ["protein_coding", "protein_coding"],
            "source_version": ["test", "test"],
        }
    )
    regions = pd.DataFrame(
        {
            "region_id": [
                "ENSG1:promoter_core",
                "RBS:cgi_1",
                "RBS:orphan_island",
                "TILE:0",
            ],
            "gene_id": ["ENSG1", None, None, None],
            "region_type": ["promoter_core", "cgi_island", "cgi_island", "cpg_tile"],
            "region_system": ["gene", "rbs", "rbs", "tbs"],
            "chromosome": ["chr1", "chr1", "chr2", "chr1"],
            "start": [90, 200, 10, 9000],
            "end": [150, 260, 80, 9100],
        }
    )
    # Loci: 0 gene, 1 near RBS→ENSG1, 2 orphan chrom RBS, 3 leftover, 4 TBS-only
    locus_index = pd.DataFrame(
        {
            "col_index": [0, 1, 2, 3, 4],
            "locus_id": [10, 20, 30, 40, 50],
            "canonical_key": [
                "GRCh38:chr1:100",
                "GRCh38:chr1:220",
                "GRCh38:chr2:40",
                "GRCh38:chr1:8000",
                "GRCh38:chr1:9050",
            ],
        }
    )
    locus_region_edges = pd.DataFrame(
        {
            "locus_id": [10, 20, 30, 50],
            "region_id": [
                "ENSG1:promoter_core",
                "RBS:cgi_1",
                "RBS:orphan_island",
                "TILE:0",
            ],
        }
    )
    n_samples = 12
    n_cols = 5
    betas = rng.uniform(0.05, 0.95, size=(n_samples, n_cols)).astype(np.float64)
    # Age correlates with gene-region column.
    ages = (40.0 + 30.0 * betas[:, 0] + rng.normal(0, 0.5, size=n_samples)).astype(np.float64)
    tissue = (betas[:, 0] > 0.5).astype(np.int64)
    sex = (betas[:, 1] > 0.5).astype(np.int64)
    sample_ids = [f"S{i}" for i in range(n_samples)]
    study_ids = np.asarray([f"ST{i // 4}" for i in range(n_samples)])
    return {
        "genes": genes,
        "regions": regions,
        "locus_index": locus_index,
        "locus_region_edges": locus_region_edges,
        "betas": betas,
        "sample_ids": sample_ids,
        "ages": ages,
        "tissue": tissue,
        "sex": sex,
        "study_ids": study_ids,
        "class_names": ["A", "B"],
    }


def _dense_cpg_features(
    beta_row: np.ndarray,
    *,
    epsilon: float = 0.001,
) -> np.ndarray:
    """M-value features [n_cols, 1] with NaN→0 (observed handled via edges)."""
    b = np.asarray(beta_row, dtype=np.float64)
    finite = np.isfinite(b)
    safe = np.where(finite, np.clip(b, epsilon, 1.0 - epsilon), 0.5)
    m = beta_to_m_value(safe, epsilon=epsilon)
    m = np.where(finite, m, 0.0)
    return m.astype(np.float32).reshape(-1, 1)


def _dense_cpg_features_batch(
    betas_rows: np.ndarray,
    *,
    col_indices: np.ndarray | None = None,
    epsilon: float = 0.001,
) -> np.ndarray:
    """M-value features [batch, n_edges, 1]; optional column subset."""
    b = np.asarray(betas_rows, dtype=np.float64)
    if b.ndim == 1:
        b = b.reshape(1, -1)
    if col_indices is not None:
        b = b[:, np.asarray(col_indices, dtype=np.int64)]
    finite = np.isfinite(b)
    safe = np.where(finite, np.clip(b, epsilon, 1.0 - epsilon), 0.5)
    m = beta_to_m_value(safe, epsilon=epsilon)
    m = np.where(finite, m, 0.0)
    return m.astype(np.float32)[..., np.newaxis]


@dataclass(frozen=True, slots=True)
class _CascadeGraphTensors:
    cpg_to_region: torch.Tensor
    region_type: torch.Tensor
    region_to_gene: torch.Tensor
    orphan_idx: torch.Tensor
    n_regions: int
    n_gene_instances: int
    edge_cols: np.ndarray


def _cascade_graph_tensors(
    assignment: CascadeAssignment,
    *,
    device: torch.device,
) -> _CascadeGraphTensors:
    cols = assignment.edge_col_index
    return _CascadeGraphTensors(
        cpg_to_region=torch.from_numpy(assignment.edge_region_index.astype(np.int64)).to(device),
        region_type=torch.from_numpy(assignment.region_type_id.astype(np.int64)).to(device),
        region_to_gene=torch.from_numpy(assignment.region_to_gene.astype(np.int64)).to(device),
        orphan_idx=torch.from_numpy(assignment.orphan_region_indices).to(device),
        n_regions=int(assignment.n_regions),
        n_gene_instances=max(int(assignment.n_genes), 1),
        edge_cols=cols,
    )


def resolve_cascade_train_batch_size(
    device: torch.device,
    *,
    n_cols: int,
    n_edges: int,
    requested: int | str | None = None,
    max_batch: int = 64,
    gpu_share: int = 1,
) -> int:
    """Resolve train/score micro-batch size: 1 on CPU; ``auto`` splits VRAM fairly.

    When several cascade jobs share one GPU (e.g. P4 ∥ P5), set ``gpu_share`` or
    env ``MBS_CASCADE_GPU_SHARE`` to the number of concurrent trainers. Reserve
    headroom for sibling jobs (encoder parity, etc.) via
    ``MBS_CASCADE_GPU_RESERVED_MIB`` (default 2048).
    """
    if requested is not None and requested != "auto":
        size = int(requested)
        if size < 1:
            raise ValueError(f"training.batch_size must be >= 1, found {size}")
        return size
    if device.type != "cuda":
        return 1
    share_env = os.environ.get("MBS_CASCADE_GPU_SHARE")
    if share_env is not None:
        gpu_share = max(1, int(share_env))
    else:
        gpu_share = max(1, int(gpu_share))
    reserved_mib = int(os.environ.get("MBS_CASCADE_GPU_RESERVED_MIB", "2048"))
    reserved_bytes = reserved_mib * 1024 * 1024
    _free_bytes, total_bytes = torch.cuda.mem_get_info(device)
    # Fair split of total VRAM — do not use full ``free`` when jobs start together.
    budget_bytes = max(
        float(total_bytes) * 0.80 - float(reserved_bytes),
        float(total_bytes) * 0.25,
    ) / float(gpu_share)
    # ponytail: linear estimate; scales with edge count (Epic/WGS) not just n_cols.
    n_edges_eff = max(int(n_edges), 1, int(n_cols))
    hidden = 64
    per_sample = float(n_edges_eff * hidden * 4 + n_cols * 4 + 48_000_000)
    batch = max(1, min(int(max_batch), int(budget_bytes / per_sample)))
    return batch


def _forward_batch(
    model: CascadeDeepSet,
    assignment: CascadeAssignment,
    betas_batch: np.ndarray,
    *,
    device: torch.device,
    graph: _CascadeGraphTensors | None = None,
) -> dict[str, torch.Tensor]:
    """Batched CpG encoder + per-sample region path; returns stacked MBS tensors."""
    betas_batch = np.asarray(betas_batch, dtype=np.float64)
    if betas_batch.ndim == 1:
        betas_batch = betas_batch.reshape(1, -1)
    batch_size = int(betas_batch.shape[0])
    graph = graph or _cascade_graph_tensors(assignment, device=device)
    cols = graph.edge_cols
    if cols.size == 0 or graph.n_regions == 0:
        mbs = torch.full(
            (batch_size, graph.n_gene_instances),
            0.5,
            device=device,
            dtype=torch.float32,
        )
        present = torch.zeros(
            batch_size, graph.n_gene_instances, dtype=torch.bool, device=device
        )
        n_orphan = int(assignment.n_orphan_rbs)
        orphan = torch.full((batch_size, n_orphan), 0.5, device=device, dtype=torch.float32)
        n_regions = int(graph.n_regions)
        rbs = torch.full((batch_size, n_regions), 0.5, device=device, dtype=torch.float32)
        rbs_present = torch.zeros(batch_size, n_regions, dtype=torch.bool, device=device)
        return {
            "mbs": mbs,
            "present": present,
            "orphan_rbs": orphan,
            "rbs": rbs,
            "rbs_present": rbs_present,
        }

    feats = _dense_cpg_features_batch(betas_batch, col_indices=cols)
    feats_t = torch.from_numpy(feats).to(device=device, dtype=torch.float32)
    n_edges = int(feats_t.shape[1])
    cpg_hidden = model.cpg_encoder(feats_t.reshape(batch_size * n_edges, 1)).view(
        batch_size, n_edges, -1
    )

    mbs_rows: list[torch.Tensor] = []
    present_rows: list[torch.Tensor] = []
    orphan_rows: list[torch.Tensor] = []
    rbs_rows: list[torch.Tensor] = []
    rbs_present_rows: list[torch.Tensor] = []
    for b in range(batch_size):
        out = model.forward_from_cpg_hidden(
            cpg_hidden[b],
            cpg_to_region=graph.cpg_to_region,
            region_type=graph.region_type,
            region_to_gene=graph.region_to_gene,
            n_regions=graph.n_regions,
            n_gene_instances=graph.n_gene_instances,
            orphan_region_indices=graph.orphan_idx,
        )
        mbs_rows.append(out["mbs"])
        present_rows.append(out["present"])
        rbs_rows.append(out["rbs"])
        rbs_present_rows.append(out["rbs_present"])
        if assignment.n_orphan_rbs:
            orphan_rows.append(out["orphan_rbs"])

    mbs = torch.stack(mbs_rows, dim=0)
    present = torch.stack(present_rows, dim=0)
    rbs = torch.stack(rbs_rows, dim=0)
    rbs_present = torch.stack(rbs_present_rows, dim=0)
    orphan = (
        torch.stack(orphan_rows, dim=0)
        if orphan_rows
        else torch.zeros(batch_size, 0, device=device, dtype=torch.float32)
    )
    if assignment.n_genes == 0:
        mbs = torch.full((batch_size, 1), 0.5, device=device, dtype=torch.float32)
        present = torch.zeros(batch_size, 1, dtype=torch.bool, device=device)
    return {
        "mbs": mbs,
        "present": present,
        "orphan_rbs": orphan,
        "rbs": rbs,
        "rbs_present": rbs_present,
    }


def _forward_sample(
    model: CascadeDeepSet,
    assignment: CascadeAssignment,
    beta_row: np.ndarray,
    *,
    device: torch.device,
) -> dict[str, torch.Tensor]:
    out = _forward_batch(model, assignment, beta_row, device=device)
    return {
        "mbs": out["mbs"][0],
        "present": out["present"][0],
        "orphan_rbs": out["orphan_rbs"][0],
    }


def score_samples(
    model: CascadeDeepSet,
    assignment: CascadeAssignment,
    betas: np.ndarray,
    *,
    device: torch.device,
    batch_size: int = 1,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return mbs, gene_present, orphan_rbs, all_rbs, all_rbs_present."""
    model.eval()
    n = betas.shape[0]
    n_genes = max(assignment.n_genes, 1)
    n_orphan = assignment.n_orphan_rbs
    n_regions = assignment.n_regions
    mbs_out = np.full((n, n_genes), 0.5, dtype=np.float32)
    present_out = np.zeros((n, n_genes), dtype=bool)
    orphan_out = np.full((n, n_orphan), 0.5, dtype=np.float32)
    rbs_out = np.full((n, n_regions), 0.5, dtype=np.float32)
    rbs_present_out = np.zeros((n, n_regions), dtype=bool)
    graph = _cascade_graph_tensors(assignment, device=device)
    step = max(1, int(batch_size))
    with torch.no_grad():
        for start in range(0, n, step):
            end = min(n, start + step)
            out = _forward_batch(
                model,
                assignment,
                betas[start:end],
                device=device,
                graph=graph,
            )
            m = out["mbs"].detach().cpu().numpy().astype(np.float32)
            p = out["present"].detach().cpu().numpy().astype(bool)
            if assignment.n_genes == 0:
                mbs_out[start:end, 0] = 0.5
                present_out[start:end, 0] = False
            else:
                mbs_out[start:end, : assignment.n_genes] = m[:, : assignment.n_genes]
                present_out[start:end, : assignment.n_genes] = p[:, : assignment.n_genes]
            if n_orphan:
                orphan_out[start:end] = (
                    out["orphan_rbs"].detach().cpu().numpy().astype(np.float32)
                )
            if n_regions:
                rbs_out[start:end] = out["rbs"].detach().cpu().numpy().astype(np.float32)
                rbs_present_out[start:end] = (
                    out["rbs_present"].detach().cpu().numpy().astype(bool)
                )
    if assignment.n_genes == 0:
        return mbs_out, present_out, orphan_out, rbs_out, rbs_present_out
    return (
        mbs_out[:, : assignment.n_genes],
        present_out[:, : assignment.n_genes],
        orphan_out,
        rbs_out,
        rbs_present_out,
    )


def _fit_direct_columns(
    *,
    betas_train: np.ndarray,
    betas_all: np.ndarray,
    direct_cols: np.ndarray,
    ages_train: np.ndarray | None,
    age_mask_train: np.ndarray | None,
    tissue_train: np.ndarray | None,
    tissue_mask_train: np.ndarray | None,
    sex_train: np.ndarray | None,
    sex_mask_train: np.ndarray | None,
    study_ids_train: np.ndarray,
    use_level1: bool = True,
) -> tuple[np.ndarray, list[str]]:
    """Fold-fit elastic-net on leftover loci; return [n_all, n_tasks] preds."""
    n_all = betas_all.shape[0]
    if direct_cols.size == 0:
        return np.zeros((n_all, 0), dtype=np.float32), []

    def _slice(betas: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        b = betas[:, direct_cols]
        obs = np.isfinite(b)
        safe = np.where(obs, b, 0.5)
        m = beta_to_m_value(safe)
        m = np.where(obs, m, np.nan)
        return m, obs

    m_tr, obs_tr = _slice(betas_train)
    m_all, obs_all = _slice(betas_all)
    z_tr, params = direct_cpg_design_matrix(m_tr, obs_tr, use_level1=use_level1)
    z_all, _ = direct_cpg_design_matrix(
        m_all, obs_all, use_level1=use_level1, level1_params=params
    )
    cols: list[np.ndarray] = []
    names: list[str] = []

    def _predict_task(y: np.ndarray, mask: np.ndarray, name: str) -> None:
        nonlocal cols, names
        if not bool(mask.any()):
            return
        try:
            fitted = fit_direct_elasticnet(
                z_tr[mask],
                obs_tr[mask],
                y[mask],
                study_ids_train[mask],
                min_studies=1,
                alpha=0.05,
            )
        except ValueError:
            return
        keep = fitted["keep_mask"]
        w = fitted["weights"]
        intercept = float(fitted["intercept"])
        x_all = np.where(obs_all, z_all, 0.0)[:, keep]
        pred = (x_all @ w + intercept).astype(np.float32)
        cols.append(pred.reshape(-1, 1))
        names.append(name)

    if ages_train is not None and age_mask_train is not None:
        _predict_task(
            np.asarray(ages_train, dtype=np.float64),
            np.asarray(age_mask_train, dtype=bool),
            "age",
        )
    if tissue_train is not None and tissue_mask_train is not None:
        _predict_task(
            np.asarray(tissue_train, dtype=np.float64),
            np.asarray(tissue_mask_train, dtype=bool),
            "tissue",
        )
    if sex_train is not None and sex_mask_train is not None:
        _predict_task(
            np.asarray(sex_train, dtype=np.float64),
            np.asarray(sex_mask_train, dtype=bool),
            "sex",
        )

    if not cols:
        return np.zeros((n_all, 0), dtype=np.float32), []
    return np.concatenate(cols, axis=1), names

def _tissue_class_weights(
    tissue_train: np.ndarray,
    tissue_mask_train: np.ndarray,
    *,
    n_classes: int,
    device: torch.device,
) -> torch.Tensor:
    """Inverse-frequency weights over classes seen in train; absent classes -> 1.0."""
    labels = tissue_train[tissue_mask_train]
    weights = np.ones(n_classes, dtype=np.float64)
    if labels.size > 0:
        counts = np.bincount(labels, minlength=n_classes).astype(np.float64)
        present = counts > 0
        if present.any():
            mean_count = counts[present].mean()
            weights[present] = mean_count / counts[present]
    return torch.tensor(weights, dtype=torch.float32, device=device)


def _evaluate_cascade_validation(
    model: CascadeDeepSet,
    heads: MultitaskHeads,
    assignment: CascadeAssignment,
    betas_val: np.ndarray,
    *,
    ages_val: np.ndarray,
    age_mask_val: np.ndarray,
    tissue_val: np.ndarray,
    tissue_mask_val: np.ndarray,
    sex_val: np.ndarray,
    sex_mask_val: np.ndarray,
    device: torch.device,
    batch_size: int = 1,
) -> dict[str, Any]:
    """Cheap proxy metrics from the model's own heads on a held-out validation slice."""
    from mbs.evaluation.metrics import (  # noqa: PLC0415
        binary_auroc_auprc,
        multiclass_metrics,
        regression_metrics,
    )

    out: dict[str, Any] = {"tissue_macro_f1": None, "age_mae": None, "sex_auroc": None}
    if betas_val.shape[0] == 0:
        return out
    mbs_v, present_v, _, _, _ = score_samples(
        model, assignment, betas_val, device=device, batch_size=batch_size
    )
    mbs_t = torch.from_numpy(mbs_v).to(device)
    present_t = torch.from_numpy(present_v).to(device)
    heads.eval()
    with torch.no_grad():
        if bool(age_mask_val.any()):
            age_pred = heads.forward_age(mbs_t, present_t).detach().cpu().numpy()
            out["age_mae"] = regression_metrics(
                ages_val[age_mask_val], age_pred[age_mask_val]
            )["mae"]
        if bool(tissue_mask_val.any()):
            tissue_logits = heads.forward_tissue(mbs_t, present_t)
            tissue_pred = tissue_logits.detach().cpu().numpy().argmax(axis=1)
            out["tissue_macro_f1"] = multiclass_metrics(
                tissue_val[tissue_mask_val], tissue_pred[tissue_mask_val]
            )["macro_f1"]
        if heads.sex_head is not None and bool(sex_mask_val.any()):
            sex_proba = (
                torch.softmax(heads.forward_sex(mbs_t, present_t), dim=-1)
                .detach()
                .cpu()
                .numpy()
            )
            yt = np.asarray(sex_val)[sex_mask_val]
            if np.unique(yt).size == 2 and sex_proba.shape[1] >= 2:
                try:
                    out["sex_auroc"] = binary_auroc_auprc(yt, sex_proba[sex_mask_val][:, 1])[
                        "auroc"
                    ]
                except ValueError:
                    out["sex_auroc"] = None
    return out


def _evaluate_mbs_e2e(
    model: CascadeDeepSet,
    heads: MultitaskHeads,
    assignment: CascadeAssignment,
    betas: np.ndarray,
    *,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    ages: np.ndarray,
    age_mask: np.ndarray,
    tissue: np.ndarray,
    tissue_mask: np.ndarray,
    sex: np.ndarray,
    sex_mask: np.ndarray,
    study_ids: np.ndarray,
    class_names: list[str],
    device: torch.device,
    batch_size: int = 1,
) -> dict[str, Any]:
    """End-to-end phenotype heads on MBS only (no late fusion); test split only."""
    test_idx_a = np.asarray(test_idx, dtype=np.int64)
    betas_te = betas[test_idx_a]
    mbs_te, present_te, _, _, _ = score_samples(
        model, assignment, betas_te, device=device, batch_size=batch_size
    )
    mbs_t = torch.from_numpy(mbs_te).to(device)
    present_t = torch.from_numpy(present_te).to(device)
    heads.eval()
    with torch.no_grad():
        age_pred = heads.forward_age(mbs_t, present_t).detach().cpu().numpy()
        tissue_logits = heads.forward_tissue(mbs_t, present_t)
        tissue_pred = tissue_logits.detach().cpu().numpy().argmax(axis=1)
        tissue_proba = torch.softmax(tissue_logits, dim=-1).detach().cpu().numpy()
        sex_pred: np.ndarray | None = None
        sex_proba: np.ndarray | None = None
        if heads.sex_head is not None:
            sex_logits = heads.forward_sex(mbs_t, present_t)
            if sex_logits is not None:
                sex_pred = sex_logits.detach().cpu().numpy().argmax(axis=1)
                sex_proba = torch.softmax(sex_logits, dim=-1).detach().cpu().numpy()
    preds: dict[str, np.ndarray | None] = {
        "age": age_pred,
        "tissue": tissue_pred,
        "tissue_proba": tissue_proba,
        "tissue_classes": np.arange(tissue_proba.shape[1], dtype=np.int64),
        "sex": sex_pred,
    }
    if sex_proba is not None:
        preds["sex_proba"] = sex_proba
        preds["sex_classes"] = np.arange(sex_proba.shape[1], dtype=np.int64)
    train_idx_a = np.asarray(train_idx, dtype=np.int64)
    tissue_valid_classes = None
    tm_tr = np.asarray(tissue_mask, dtype=bool)[train_idx_a]
    if tm_tr.any():
        tissue_valid_classes = set(
            np.asarray(tissue, dtype=np.int64)[train_idx_a][tm_tr].tolist()
        )
    metrics = evaluate_multitask_predictions(
        preds=preds,
        age=ages[test_idx_a],
        age_mask=age_mask[test_idx_a],
        tissue=tissue[test_idx_a],
        tissue_mask=tissue_mask[test_idx_a],
        sex=sex[test_idx_a],
        sex_mask=sex_mask[test_idx_a],
        study_ids=study_ids[test_idx_a],
        tissue_class_names=list(class_names) if class_names else None,
        tissue_valid_classes=tissue_valid_classes,
    )
    return {
        "metrics": metrics,
        "evaluation": "mbs_e2e",
        "eval_split": "test",
        "n_eval_samples": int(test_idx_a.size),
        "n_score_features": int(mbs_te.shape[1]),
    }


def _evaluate_fusion_mode(
    blocks: dict[str, np.ndarray],
    *,
    mode: FusionBlockMode,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    ages: np.ndarray,
    age_mask: np.ndarray,
    tissue: np.ndarray,
    tissue_mask: np.ndarray,
    sex: np.ndarray,
    sex_mask: np.ndarray,
    study_ids: np.ndarray,
    class_names: list[str],
    fusion: dict[str, Any] | None,
) -> dict[str, Any]:
    x = fusion_feature_matrix(blocks, mode=mode)
    x_tr = x[train_idx]
    x_te = x[test_idx]
    out = evaluate_late_fusion(
        scores_train=x_tr,
        scores_test=x_te,
        age_train=ages[train_idx],
        age_mask_train=age_mask[train_idx],
        tissue_train=tissue[train_idx],
        tissue_mask_train=tissue_mask[train_idx],
        sex_train=sex[train_idx],
        sex_mask_train=sex_mask[train_idx],
        age_test=ages[test_idx],
        age_mask_test=age_mask[test_idx],
        tissue_test=tissue[test_idx],
        tissue_mask_test=tissue_mask[test_idx],
        sex_test=sex[test_idx],
        sex_mask_test=sex_mask[test_idx],
        study_ids_test=study_ids[test_idx],
        tissue_class_names=list(class_names) if class_names else None,
        fusion=fusion,
    )
    out["evaluation"] = f"fusion_{mode}"
    out["fusion_block_mode"] = mode
    return out


def _evaluate_mbs_enet(
    blocks: dict[str, np.ndarray],
    *,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    ages: np.ndarray,
    age_mask: np.ndarray,
    tissue: np.ndarray,
    tissue_mask: np.ndarray,
    sex: np.ndarray,
    sex_mask: np.ndarray,
    study_ids: np.ndarray,
    class_names: list[str],
) -> dict[str, Any]:
    """Elastic-net readout on frozen MBS columns only (same encoder as e2e / linear probe)."""
    x = fusion_feature_matrix(blocks, mode="mbs_only")
    out = run_elasticnet_multitask(
        x_train=x[train_idx],
        x_test=x[test_idx],
        age_train=ages[train_idx],
        age_mask_train=age_mask[train_idx],
        tissue_train=tissue[train_idx],
        tissue_mask_train=tissue_mask[train_idx],
        sex_train=sex[train_idx],
        sex_mask_train=sex_mask[train_idx],
        age_test=ages[test_idx],
        age_mask_test=age_mask[test_idx],
        tissue_test=tissue[test_idx],
        tissue_mask_test=tissue_mask[test_idx],
        sex_test=sex[test_idx],
        sex_mask_test=sex_mask[test_idx],
        study_ids_test=study_ids[test_idx],
        tissue_class_names=list(class_names) if class_names else None,
    )
    out["evaluation"] = "mbs_enet"
    out["eval_split"] = "test"
    out["n_eval_samples"] = int(np.asarray(test_idx).size)
    out["n_score_features"] = int(out.get("n_features", x.shape[1]))
    return out


def _evaluations_incomplete(metrics: dict[str, Any]) -> bool:
    ev = metrics.get("evaluations")
    if not isinstance(ev, dict):
        return True
    if "mbs_e2e" not in ev or "fusion_full" not in ev:
        return True
    e2e = ev.get("mbs_e2e")
    if not isinstance(e2e, dict) or e2e.get("eval_split") != "test":
        return True
    # Gene-linked Stage A: need frozen RBS linear probe. Elastic-net (mbs_enet /
    # rbs_enet) is post-hoc and must not block the GPU screen queue.
    if bool(metrics.get("gene_linked_only")) and "rbs_linear_probe" not in ev:
        return True
    return False


def train_cascade_on_arrays(
    *,
    assignment: CascadeAssignment,
    betas: np.ndarray,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    ages: np.ndarray,
    tissue: np.ndarray,
    sex: np.ndarray,
    study_ids: np.ndarray,
    sample_ids: list[str],
    class_names: list[str],
    out_dir: Path,
    max_epochs: int = 8,
    seed: int = 42,
    device_str: str = "cpu",
    lr: float = 1e-2,
    age_mask: np.ndarray | None = None,
    tissue_mask: np.ndarray | None = None,
    sex_mask: np.ndarray | None = None,
    age_seed_mask: np.ndarray | torch.Tensor | None = None,
    tissue_seed_mask: np.ndarray | torch.Tensor | None = None,
    sex_seed_mask: np.ndarray | torch.Tensor | None = None,
    checkpoint_selection_mode: str = TISSUE_PRIMARY_SELECTION,
    cpg_hidden_dim: int = 64,
    region_hidden_dim: int = 32,
    dropout: float = 0.1,
    cpg_pool: PoolName = "max",
    region_pool: PoolName = "max",
    gene_aggregation: str = "scalar_rbs",
    gene_allocation_policy: str = "legacy_nearest",
    skip_if_done: bool = False,
    val_idx: np.ndarray | None = None,
    age_loss_weight: float = 1.0,
    tissue_loss_weight: float = 1.0,
    sex_loss_weight: float = 1.0,
    early_stopping_patience: int | None = None,
    early_stopping_min_delta: float = 0.0,
    fusion: dict[str, Any] | None = None,
    gene_linked_only: bool = False,
    primary_evaluation: PrimaryEvaluation = "late_fusion",
    extra_fusion_modes: tuple[FusionBlockMode, ...] = (),
    locus_ids: list[str] | None = None,
    eval_only: bool = False,
    train_batch_size: int | str | None = "auto",
    gpu_share: int = 1,
    include_mbs_enet: bool = True,
) -> dict[str, Any]:
    """Train CascadeDeepSet + MBS heads; write scores; evaluate; return metrics.

    ``include_mbs_enet``: when False, skip inline ``mbs_enet`` / ``rbs_enet`` (Stage A
    screen default — run ``scripts/eval_mbs_enet_from_scores.py`` post-hoc).
    """
    score_dir = out_dir / "scores"
    manifest_path = score_dir / "score_manifest.json"
    metrics_path = out_dir / "metrics.json"
    prior_metrics: dict[str, Any] | None = None
    if metrics_path.is_file():
        prior_metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
    if (
        skip_if_done
        and manifest_path.is_file()
        and metrics_path.is_file()
        and prior_metrics is not None
        and not _evaluations_incomplete(prior_metrics)
    ):
        cached = dict(prior_metrics or {})
        cached["skipped"] = True
        cached["score_dir"] = str(score_dir)
        return cached

    if early_stopping_patience is not None and early_stopping_patience < 1:
        raise ValueError("early_stopping_patience must be >= 1 when enabled")
    if early_stopping_min_delta < 0:
        raise ValueError("early_stopping_min_delta must be non-negative")
    if primary_evaluation not in ("late_fusion", "mbs_e2e"):
        raise ValueError(f"unsupported primary_evaluation: {primary_evaluation!r}")
    if checkpoint_selection_mode not in (AGE_PRIMARY_SELECTION, TISSUE_PRIMARY_SELECTION):
        raise ValueError(f"unsupported checkpoint_selection_mode: {checkpoint_selection_mode!r}")

    if gene_linked_only:
        assignment = assignment_gene_linked_only(assignment)
    gene_cols = gene_linked_col_index(assignment)

    _set_seed(seed)
    device = resolve_device(device_str, require_cuda=False)
    n_region_types = max(len(assignment.region_types), 1)
    n_genes = max(assignment.n_genes, 1)
    model = CascadeDeepSet(
        1,
        n_region_types,
        cpg_hidden_dim=int(cpg_hidden_dim),
        region_hidden_dim=int(region_hidden_dim),
        cpg_pool=cpg_pool,
        region_pool=region_pool,
        gene_aggregation=gene_aggregation,  # type: ignore[arg-type]
        dropout=float(dropout),
        activation="gelu",
        layer_norm=True,
    )
    age_seed_t = _seed_mask_tensor(age_seed_mask)
    tissue_seed_t = _seed_mask_tensor(tissue_seed_mask)
    sex_seed_t = _seed_mask_tensor(sex_seed_mask)
    # Fails closed (ValueError) before training if any provided mask is undersized.
    heads = MultitaskHeads(
        n_genes,
        max(len(class_names), 2),
        sex_enabled=True,
        age_seed_mask=age_seed_t,
        tissue_seed_mask=tissue_seed_t,
        sex_seed_mask=sex_seed_t,
    )
    seed_mask_meta = {
        trait: meta
        for trait, meta in (
            ("age", _seed_mask_meta(age_seed_t)),
            ("tissue", _seed_mask_meta(tissue_seed_t)),
            ("sex", _seed_mask_meta(sex_seed_t)),
        )
        if meta is not None
    }
    model.to(device)
    heads.to(device)
    opt = torch.optim.Adam(list(model.parameters()) + list(heads.parameters()), lr=lr)

    train_idx = np.asarray(train_idx, dtype=np.int64)
    test_idx = np.asarray(test_idx, dtype=np.int64)
    n = len(sample_ids)
    age_mask_a = (
        np.ones(n, dtype=bool) if age_mask is None else np.asarray(age_mask, dtype=bool)
    )
    tissue_mask_a = (
        np.ones(n, dtype=bool)
        if tissue_mask is None
        else np.asarray(tissue_mask, dtype=bool)
    )
    sex_mask_a = (
        np.ones(n, dtype=bool) if sex_mask is None else np.asarray(sex_mask, dtype=bool)
    )
    tissue_class_weights = _tissue_class_weights(
        tissue[train_idx],
        tissue_mask_a[train_idx],
        n_classes=max(len(class_names), 2),
        device=device,
    )
    batch_size = resolve_cascade_train_batch_size(
        device,
        n_cols=int(betas.shape[1]),
        n_edges=int(assignment.edge_col_index.size),
        requested=train_batch_size,
        gpu_share=gpu_share,
    )
    graph_tensors = _cascade_graph_tensors(assignment, device=device)
    share_note = ""
    if device.type == "cuda":
        share_env = os.environ.get("MBS_CASCADE_GPU_SHARE")
        effective_share = max(1, int(share_env)) if share_env else max(1, int(gpu_share))
        if effective_share > 1:
            share_note = f" gpu_share={effective_share}"
    print(
        f"[cascade] {out_dir.name} train_batch_size={batch_size} "
        f"n_edges={assignment.edge_col_index.size} device={device.type}{share_note}",
        flush=True,
    )

    val_idx_a = None if val_idx is None else np.asarray(val_idx, dtype=np.int64)
    has_val = val_idx_a is not None and val_idx_a.size > 0
    out_dir.mkdir(parents=True, exist_ok=True)
    ckpt_path = out_dir / "best.pt"
    # Re-score / re-eval without retraining when a checkpoint already exists but
    # Stage A diagnostics (e.g. rbs_linear_probe) are missing.
    if (
        not eval_only
        and ckpt_path.is_file()
        and prior_metrics is not None
        and _evaluations_incomplete(prior_metrics)
    ):
        eval_only = True
        print(
            f"[cascade] {out_dir.name} eval_only auto: checkpoint present, "
            "refreshing incomplete evaluations",
            flush=True,
        )
    if eval_only and not ckpt_path.is_file():
        raise FileNotFoundError(f"eval_only requires checkpoint: {ckpt_path}")

    def _save_checkpoint() -> None:
        torch.save(
            {
                "model": model.state_dict(),
                "heads": heads.state_dict(),
                "seed": seed,
                "max_epochs": max_epochs,
                "cpg_hidden_dim": cpg_hidden_dim,
                "region_hidden_dim": region_hidden_dim,
                "cpg_pool": cpg_pool,
                "region_pool": region_pool,
                "gene_aggregation": gene_aggregation,
            },
            ckpt_path,
        )

    age_primary = checkpoint_selection_mode == AGE_PRIMARY_SELECTION
    rank_fn = _validation_rank_age_primary if age_primary else _validation_rank
    early_stop_monitor = "validation_age_mae" if age_primary else "validation_tissue_macro_f1"

    best_rank: tuple[float, ...] | None = None
    best_epoch = -1
    val_history: list[dict[str, Any]] = []
    best_early_stop_f1 = -float("inf")
    epochs_without_tissue_improvement = 0
    epochs_completed = 0
    stopped_early = False
    stop_epoch: int | None = None

    if eval_only:
        ckpt = torch.load(ckpt_path, map_location=device)
        model.load_state_dict(ckpt["model"])
        heads.load_state_dict(ckpt["heads"])
        if prior_metrics and isinstance(prior_metrics.get("checkpoint_selection"), dict):
            checkpoint_selection = dict(prior_metrics["checkpoint_selection"])
        else:
            checkpoint_selection = {"selection": "eval_only_reload", "best_epoch": ckpt.get("epoch")}
    else:
        for _epoch in range(max_epochs):
            epochs_completed = _epoch + 1
            model.train()
            heads.train()
            order = train_idx.copy()
            np.random.shuffle(order)
            batch_starts = range(0, int(order.size), batch_size)
            log_this_epoch = (_epoch + 1) % max(1, max_epochs // 5) == 0 or _epoch == 0
            logged_grad_norms = False
            for start in batch_starts:
                batch_idx = order[start : start + batch_size]
                active = [
                    int(i)
                    for i in batch_idx.tolist()
                    if bool(age_mask_a[i]) or bool(tissue_mask_a[i]) or bool(sex_mask_a[i])
                ]
                if not active:
                    continue
                active_a = np.asarray(active, dtype=np.int64)
                out = _forward_batch(
                    model,
                    assignment,
                    betas[active_a],
                    device=device,
                    graph=graph_tensors,
                )
                mbs = out["mbs"]
                present = out["present"]
                if assignment.n_genes == 0:
                    mbs = torch.full((len(active), 1), 0.5, device=device)
                    present = torch.zeros(len(active), 1, dtype=torch.bool, device=device)
                loss = torch.zeros((), device=device)
                age_term: torch.Tensor | None = None
                tissue_term: torch.Tensor | None = None
                sex_term: torch.Tensor | None = None
                age_m = age_mask_a[active_a]
                if bool(age_m.any()):
                    age_t = torch.tensor(ages[active_a][age_m], device=device, dtype=torch.float32)
                    age_pred = heads.forward_age(mbs[age_m], present[age_m])
                    age_term = F.huber_loss(age_pred, age_t)
                    loss = loss + age_loss_weight * age_term
                tissue_m = tissue_mask_a[active_a]
                if bool(tissue_m.any()):
                    tissue_t = torch.tensor(
                        tissue[active_a][tissue_m], device=device, dtype=torch.long
                    )
                    tissue_pred = heads.forward_tissue(mbs[tissue_m], present[tissue_m])
                    tissue_term = F.cross_entropy(
                        tissue_pred, tissue_t, weight=tissue_class_weights
                    )
                    loss = loss + tissue_loss_weight * tissue_term
                sex_m = sex_mask_a[active_a]
                if bool(sex_m.any()):
                    sex_t = torch.tensor(sex[active_a][sex_m], device=device, dtype=torch.long)
                    sex_pred = heads.forward_sex(mbs[sex_m], present[sex_m])
                    if sex_pred is not None:
                        sex_term = F.cross_entropy(sex_pred, sex_t)
                        loss = loss + sex_loss_weight * sex_term
                # Age-primary screen: log unweighted trait losses + encoder ||grad||
                # per trait once per logged epoch (first labelled batch).
                if log_this_epoch and not logged_grad_norms:
                    enc_norms: dict[str, float] = {}
                    for trait_name, term in (
                        ("age", age_term),
                        ("tissue", tissue_term),
                        ("sex", sex_term),
                    ):
                        if term is None:
                            continue
                        opt.zero_grad(set_to_none=True)
                        term.backward(retain_graph=True)
                        enc_norms[trait_name] = _encoder_grad_norm(model)
                    uw = {
                        "age_loss": None if age_term is None else float(age_term.detach()),
                        "tissue_loss": None if tissue_term is None else float(tissue_term.detach()),
                        "sex_loss": None if sex_term is None else float(sex_term.detach()),
                    }
                    print(
                        f"[cascade] {out_dir.name} epoch {_epoch + 1}/{max_epochs} "
                        f"unweighted={uw} encoder_grad_norm={enc_norms}",
                        flush=True,
                    )
                    logged_grad_norms = True
                opt.zero_grad(set_to_none=True)
                loss.backward()
                opt.step()
            if log_this_epoch and not logged_grad_norms:
                print(
                    f"[cascade] {out_dir.name} epoch {_epoch + 1}/{max_epochs}",
                    flush=True,
                )
            if has_val and val_idx_a is not None:
                val_metrics = _evaluate_cascade_validation(
                    model,
                    heads,
                    assignment,
                    betas[val_idx_a],
                    ages_val=ages[val_idx_a],
                    age_mask_val=age_mask_a[val_idx_a],
                    tissue_val=tissue[val_idx_a],
                    tissue_mask_val=tissue_mask_a[val_idx_a],
                    sex_val=sex[val_idx_a],
                    sex_mask_val=sex_mask_a[val_idx_a],
                    device=device,
                    batch_size=batch_size,
                )
                rank = rank_fn(val_metrics)
                val_history.append({"epoch": _epoch + 1, "rank": list(rank), **val_metrics})
                if best_rank is None or rank > best_rank:
                    best_rank = rank
                    best_epoch = _epoch + 1
                    _save_checkpoint()

                # Monitor the primary metric as "higher is better": tissue macro-F1
                # directly, or negated age MAE under age-primary selection.
                if age_primary:
                    monitor_raw = val_metrics.get("age_mae")
                    monitor_val = None if monitor_raw is None else -float(monitor_raw)
                else:
                    monitor_raw = val_metrics.get("tissue_macro_f1")
                    monitor_val = None if monitor_raw is None else float(monitor_raw)
                if early_stopping_patience is not None and monitor_val is not None:
                    if monitor_val > best_early_stop_f1 + early_stopping_min_delta:
                        best_early_stop_f1 = monitor_val
                        epochs_without_tissue_improvement = 0
                    else:
                        epochs_without_tissue_improvement += 1
                    if epochs_without_tissue_improvement >= early_stopping_patience:
                        stopped_early = True
                        stop_epoch = _epoch + 1
            if stopped_early:
                print(
                    f"[cascade] {out_dir.name} early stop at epoch {stop_epoch}; "
                    f"best validation {early_stop_monitor}={best_early_stop_f1:.4f}",
                    flush=True,
                )
                break

        checkpoint_selection: dict[str, Any] = {
            "has_validation": has_val,
            "n_val": int(val_idx_a.size) if val_idx_a is not None else 0,
            "max_epochs": max_epochs,
            "epochs_completed": epochs_completed,
            "val_history": val_history,
            "early_stopping": {
                "enabled": early_stopping_patience is not None,
                "monitor": early_stop_monitor,
                "patience": early_stopping_patience,
                "min_delta": early_stopping_min_delta,
                "stopped_early": stopped_early,
                "stop_epoch": stop_epoch,
            },
        }
        if has_val and best_epoch > 0:
            checkpoint_selection["best_epoch"] = best_epoch
            checkpoint_selection["selection"] = checkpoint_selection_mode
            checkpoint_selection["best_validation_metrics"] = next(
                row for row in val_history if int(row["epoch"]) == best_epoch
            )
            # Reload the best-validation checkpoint (may not be the final epoch).
            ckpt = torch.load(ckpt_path, map_location=device)
            model.load_state_dict(ckpt["model"])
            heads.load_state_dict(ckpt["heads"])
        else:
            checkpoint_selection["best_epoch"] = max_epochs
            checkpoint_selection["selection"] = "final_epoch_no_validation"
            _save_checkpoint()

    if eval_only and manifest_path.is_file():
        # Re-attach diagnostic gene-linked RBS if an older score dir lacks it.
        if not (score_dir / "all_gene_rbs.zarr").exists():
            mbs_all, present_all, orphan_all, rbs_all, rbs_present_all = score_samples(
                model, assignment, betas, device=device, batch_size=batch_size
            )
            gene_linked_region_mask = assignment.region_to_gene >= 0
            gene_linked_region_indices = np.flatnonzero(gene_linked_region_mask).astype(
                np.int64
            )
            all_gene_rbs = (
                rbs_all[:, gene_linked_region_indices]
                if gene_linked_region_indices.size
                else np.zeros((rbs_all.shape[0], 0), dtype=np.float32)
            )
            all_gene_rbs_present = (
                rbs_present_all[:, gene_linked_region_indices]
                if gene_linked_region_indices.size
                else np.zeros((rbs_present_all.shape[0], 0), dtype=bool)
            )
            _write_array(score_dir / "all_gene_rbs.zarr", all_gene_rbs)
            _write_array(
                score_dir / "all_gene_rbs_present.zarr",
                all_gene_rbs_present.astype(np.uint8),
            )
            pd.DataFrame(
                {
                    "region_id": [
                        assignment.region_ids[int(i)]
                        for i in gene_linked_region_indices.tolist()
                    ],
                    "gene_id": [
                        assignment.gene_ids[int(assignment.region_to_gene[int(i)])]
                        if int(assignment.region_to_gene[int(i)]) >= 0
                        else None
                        for i in gene_linked_region_indices.tolist()
                    ],
                    "region_type": [
                        assignment.region_types[int(assignment.region_type_id[int(i)])]
                        if int(assignment.region_type_id[int(i)])
                        < len(assignment.region_types)
                        else "unknown"
                        for i in gene_linked_region_indices.tolist()
                    ],
                    "column_index": np.arange(
                        gene_linked_region_indices.size, dtype=np.int64
                    ),
                    "allocation_policy": str(gene_allocation_policy),
                }
            ).to_parquet(score_dir / "all_gene_region_index.parquet", index=False)
            del mbs_all, present_all, orphan_all
    else:
        mbs_all, present_all, orphan_all, rbs_all, rbs_present_all = score_samples(
            model, assignment, betas, device=device, batch_size=batch_size
        )
        direct_all, direct_names = _fit_direct_columns(
            betas_train=betas[train_idx],
            betas_all=betas,
            direct_cols=assignment.direct_col_index,
            ages_train=ages[train_idx],
            age_mask_train=age_mask_a[train_idx],
            tissue_train=tissue[train_idx],
            tissue_mask_train=tissue_mask_a[train_idx],
            sex_train=sex[train_idx],
            sex_mask_train=sex_mask_a[train_idx],
            study_ids_train=study_ids[train_idx],
            use_level1=True,
        )

        gene_ids = list(assignment.gene_ids) if assignment.gene_ids else ["__none__"]
        orphan_ids = list(assignment.orphan_region_ids)
        gene_linked_region_mask = assignment.region_to_gene >= 0
        gene_linked_region_indices = np.flatnonzero(gene_linked_region_mask).astype(np.int64)
        all_gene_rbs = (
            rbs_all[:, gene_linked_region_indices]
            if gene_linked_region_indices.size
            else np.zeros((rbs_all.shape[0], 0), dtype=np.float32)
        )
        all_gene_rbs_present = (
            rbs_present_all[:, gene_linked_region_indices]
            if gene_linked_region_indices.size
            else np.zeros((rbs_present_all.shape[0], 0), dtype=bool)
        )
        all_gene_region_ids = [assignment.region_ids[int(i)] for i in gene_linked_region_indices.tolist()]
        all_gene_region_gene_ids = [
            assignment.gene_ids[int(assignment.region_to_gene[int(i)])]
            if int(assignment.region_to_gene[int(i)]) >= 0
            else None
            for i in gene_linked_region_indices.tolist()
        ]
        all_gene_region_types = [
            assignment.region_types[int(assignment.region_type_id[int(i)])]
            if int(assignment.region_type_id[int(i)]) < len(assignment.region_types)
            else "unknown"
            for i in gene_linked_region_indices.tolist()
        ]
        direct_cpg_all: np.ndarray | None = None
        direct_locus_ids: list[str] | None = None
        if assignment.direct_col_index.size and locus_ids is not None:
            dcols = assignment.direct_col_index
            block = betas[:, dcols]
            obs = np.isfinite(block)
            safe = np.where(obs, block, 0.5)
            m_vals = beta_to_m_value(safe)
            direct_cpg_all = np.where(obs, m_vals, np.nan).astype(np.float32)
            direct_locus_ids = [str(locus_ids[int(c)]) for c in dcols.tolist()]
        write_cascade_score_dir(
            score_dir,
            sample_ids=sample_ids,
            gene_ids=gene_ids,
            orphan_region_ids=orphan_ids,
            mbs=mbs_all if assignment.n_genes else mbs_all[:, :1],
            gene_present=present_all if assignment.n_genes else present_all[:, :1],
            orphan_rbs=orphan_all,
            direct_contrib=direct_all,
            direct_task_names=direct_names,
            fold_id=str(out_dir.name),
            restart_id=str(seed),
            direct_cpg=direct_cpg_all,
            direct_locus_ids=direct_locus_ids,
            all_gene_rbs=all_gene_rbs,
            all_gene_rbs_present=all_gene_rbs_present,
            all_gene_region_ids=all_gene_region_ids,
            all_gene_region_gene_ids=all_gene_region_gene_ids,
            all_gene_region_types=all_gene_region_types,
            allocation_policy=str(gene_allocation_policy),
            extra_manifest={
                "cpg_pool": cpg_pool,
                "region_pool": region_pool,
                "gene_aggregation": gene_aggregation,
                "checkpoint_selection_mode": checkpoint_selection_mode,
                **({"seed_masks": seed_mask_meta} if seed_mask_meta else {}),
            },
        )

    blocks = load_cascade_score_blocks(score_dir)
    if "tbs" in blocks:
        raise RuntimeError("TBS leaked into 7F fusion")

    evaluations: dict[str, Any] = {}
    evaluations["mbs_e2e"] = _evaluate_mbs_e2e(
        model,
        heads,
        assignment,
        betas,
        train_idx=train_idx,
        test_idx=test_idx,
        ages=ages,
        age_mask=age_mask_a,
        tissue=tissue,
        tissue_mask=tissue_mask_a,
        sex=sex,
        sex_mask=sex_mask_a,
        study_ids=study_ids,
        class_names=list(class_names),
        device=device,
        batch_size=batch_size,
    )
    evaluations["mbs_linear_probe"] = _evaluate_fusion_mode(
        blocks,
        mode="mbs_only",
        train_idx=train_idx,
        test_idx=test_idx,
        ages=ages,
        age_mask=age_mask_a,
        tissue=tissue,
        tissue_mask=tissue_mask_a,
        sex=sex,
        sex_mask=sex_mask_a,
        study_ids=study_ids,
        class_names=list(class_names),
        fusion=fusion,
    )
    if include_mbs_enet:
        evaluations["mbs_enet"] = _evaluate_mbs_enet(
            blocks,
            train_idx=train_idx,
            test_idx=test_idx,
            ages=ages,
            age_mask=age_mask_a,
            tissue=tissue,
            tissue_mask=tissue_mask_a,
            sex=sex,
            sex_mask=sex_mask_a,
            study_ids=study_ids,
            class_names=list(class_names),
        )
    else:
        print("[cascade] skipping inline mbs_enet (post-hoc)", flush=True)
    if "all_gene_rbs" in blocks and blocks["all_gene_rbs"].shape[1] > 0:
        rbs_blocks = {"mbs": blocks["all_gene_rbs"]}
        evaluations["rbs_linear_probe"] = _evaluate_fusion_mode(
            rbs_blocks,
            mode="mbs_only",
            train_idx=train_idx,
            test_idx=test_idx,
            ages=ages,
            age_mask=age_mask_a,
            tissue=tissue,
            tissue_mask=tissue_mask_a,
            sex=sex,
            sex_mask=sex_mask_a,
            study_ids=study_ids,
            class_names=list(class_names),
            fusion=fusion,
        )
        evaluations["rbs_linear_probe"]["evaluation"] = "rbs_linear_probe"
        if include_mbs_enet:
            evaluations["rbs_enet"] = _evaluate_mbs_enet(
                rbs_blocks,
                train_idx=train_idx,
                test_idx=test_idx,
                ages=ages,
                age_mask=age_mask_a,
                tissue=tissue,
                tissue_mask=tissue_mask_a,
                sex=sex,
                sex_mask=sex_mask_a,
                study_ids=study_ids,
                class_names=list(class_names),
            )
            evaluations["rbs_enet"]["evaluation"] = "rbs_enet"

    fusion_modes: list[FusionBlockMode] = ["full"]
    for mode in extra_fusion_modes:
        if mode not in fusion_modes:
            fusion_modes.append(mode)
    for mode in fusion_modes:
        key = f"fusion_{mode}"
        evaluations[key] = _evaluate_fusion_mode(
            blocks,
            mode=mode,
            train_idx=train_idx,
            test_idx=test_idx,
            ages=ages,
            age_mask=age_mask_a,
            tissue=tissue,
            tissue_mask=tissue_mask_a,
            sex=sex,
            sex_mask=sex_mask_a,
            study_ids=study_ids,
            class_names=list(class_names),
            fusion=fusion,
        )

    primary_key = "mbs_e2e" if primary_evaluation == "mbs_e2e" else "fusion_full"
    primary_blob = evaluations[primary_key]
    fused: dict[str, Any] = {
        "metrics": primary_blob["metrics"],
        "primary_evaluation": primary_evaluation,
        "evaluations": evaluations,
        "gene_linked_only": bool(gene_linked_only),
        "include_mbs_enet": bool(include_mbs_enet),
        "n_gene_cols": int(gene_cols.size),
        "n_orphan_rbs": int(assignment.n_orphan_rbs),
        "n_direct": int(assignment.n_direct),
        "n_genes": int(assignment.n_genes),
        "tbs_arm": False,
        "score_dir": str(score_dir),
        "fusion_n_features": int(
            evaluations.get("fusion_full", {}).get("n_score_features", 0)
        ),
        "skipped": False,
        "checkpoint": str(ckpt_path),
        "checkpoint_selection": checkpoint_selection,
        "pooling": {"cpg_to_region": cpg_pool, "region_to_gene": region_pool},
        "gene_aggregation": gene_aggregation,
        "gene_allocation": gene_allocation_policy,
        "checkpoint_selection_mode": checkpoint_selection_mode,
    }
    if seed_mask_meta:
        fused["seed_masks"] = seed_mask_meta
    write_json(metrics_path, fused)
    return fused


def run_cascade_fixture(
    *,
    project_root: Path,
    artifact_root: Path,
    run_id: str = "stage0-7f-fixture",
    seed: int = 42,
    max_epochs: int = 12,
    device_str: str = "cpu",
) -> CascadeTrainResult:
    """End-to-end fixture cascade + inspection report."""
    tables = make_synthetic_cascade_tables(seed=seed)
    assignment = build_cascade_assignment(
        locus_index=tables["locus_index"],
        locus_region_edges=tables["locus_region_edges"],
        regions=tables["regions"],
        genes=tables["genes"],
    )
    n = len(tables["sample_ids"])
    train_idx = np.arange(0, max(2, (n * 2) // 3), dtype=np.int64)
    test_idx = np.arange(train_idx[-1] + 1, n, dtype=np.int64)
    if test_idx.size == 0:
        test_idx = train_idx.copy()

    run_root = run_dir(artifact_root, run_id)
    run_root.mkdir(parents=True, exist_ok=True)
    metrics = train_cascade_on_arrays(
        assignment=assignment,
        betas=tables["betas"],
        train_idx=train_idx,
        test_idx=test_idx,
        ages=tables["ages"],
        tissue=tables["tissue"],
        sex=tables["sex"],
        study_ids=tables["study_ids"],
        sample_ids=tables["sample_ids"],
        class_names=tables["class_names"],
        out_dir=run_root,
        max_epochs=max_epochs,
        seed=seed,
        device_str=device_str,
        cpg_hidden_dim=16,
        region_hidden_dim=8,
        dropout=0.0,
    )
    report_dir = project_root / "reports" / "inspection" / "stage0_7f_rbs_gene_direct"
    report_dir.mkdir(parents=True, exist_ok=True)
    summary = {
        "milestone": "7F",
        "topology": "rbs_gene_direct",
        "tbs_arm": False,
        "split_id": "fixture",
        "assignment": {
            "n_genes": assignment.n_genes,
            "n_regions": assignment.n_regions,
            "n_orphan_rbs": assignment.n_orphan_rbs,
            "n_direct": assignment.n_direct,
            "orphan_region_ids": assignment.orphan_region_ids,
            "allocated_gene_id": assignment.allocated_gene_id,
            "direct_col_index": assignment.direct_col_index.tolist(),
        },
        "metrics": metrics,
        "score_dir": metrics.get("score_dir"),
        "adr": "0009",
    }
    write_json(report_dir / "summary.json", summary)
    analysis = f"""# Stage 0 Milestone 7F — RBS→gene + direct leftover

Status: **done** for topology acceptance (assignment + trainer + saved-score
fusion + inspection). Full-budget methylation-only bake-off is **7G**.

## Topology (ADR 0009)

```text
CpG → typed region (gene roles | RBS) → RBS score
        ├─ allocated to gene (typed and/or nearest-gene) → MBS
        └─ no gene allocation → orphan RBS
CpG with no typed region (incl. former TBS-only) → direct
late fusion: [orphan RBS | MBS | direct] → linear heads
```

- **No TBS arm** in the model matrix or score export.
- Nearest-gene applies only to typed **RBS** allocation, never leftover CpGs.

## Assignment

Nearest-gene RBS allocation: `{assignment.allocated_gene_id}`
Direct columns (not nearest-gene as CpGs): `{assignment.direct_col_index.tolist()}`

## Fusion

Features = `[orphan_rbs | mbs | direct]` (n={metrics.get("fusion_n_features")}).
`tbs_arm={metrics.get("tbs_arm")}`.

Scores: `{metrics.get("score_dir")}`.

## Metrics (fixture holdout)

```json
{json.dumps(metrics.get("metrics", metrics), indent=2, default=str)}
```

Hub re-run on frozen `hub-ats-7e-3fold-v1` uses the same cascade path
(`mbs train cascade` without `--overfit-fixture`).
"""
    (report_dir / "analysis.md").write_text(analysis, encoding="utf-8")
    return CascadeTrainResult(
        metrics=summary, score_dir=Path(str(metrics["score_dir"])), report_dir=report_dir
    )


def run_cascade_hub(
    *,
    project_root: Path,
    data_root: Path,
    artifact_root: Path,
    config: dict[str, Any],
    run_id: str,
    device_str: str = "cpu",
    max_folds: int | None = None,
    max_train_samples: int | None = None,
    report_dir: Path | None = None,
    skip_if_done: bool = True,
    eval_only: bool = False,
) -> CascadeTrainResult:
    """Train cascade on frozen 7E folds; write scores + report."""
    pilot = config.get("pilot", {})
    matrix_id = str(pilot.get("matrix_id", "matrix-hub-age-tissue-sex-full-v1"))
    graph_id = str(pilot.get("graph_id", "graph-grch38-gencode38-cgi-tile-v2"))
    split_id = str(config.get("split_id", DEFAULT_SPLIT_ID))
    cv_budget = config.get("cv_budget", {})
    max_loci = int(cv_budget.get("max_loci", config.get("training", {}).get("max_loci", 8192)))
    max_epochs = int(cv_budget.get("max_epochs", config.get("training", {}).get("max_epochs", 2)))
    seed = int(config.get("experiment", {}).get("seed", 42))
    model_cfg = config.get("model", {})
    cpg_hidden = int(model_cfg.get("cpg_hidden_dim", 64))
    region_hidden = int(model_cfg.get("region_hidden_dim", 32))
    dropout = float(model_cfg.get("dropout", 0.1))
    pool_names = {"sum", "mean", "sqrt_sum", "max"}
    default_pool = str(model_cfg.get("pooling", "max"))
    cpg_pool_name = str(model_cfg.get("cpg_pool", default_pool))
    region_pool_name = str(model_cfg.get("region_pool", default_pool))
    if cpg_pool_name not in pool_names or region_pool_name not in pool_names:
        raise ValueError(
            "model pooling must be one of sum, mean, sqrt_sum, max; "
            f"found cpg_pool={cpg_pool_name!r}, region_pool={region_pool_name!r}"
        )
    cpg_pool = cast(PoolName, cpg_pool_name)
    region_pool = cast(PoolName, region_pool_name)
    gene_aggregation = str(model_cfg.get("gene_aggregation", "scalar_rbs"))
    if gene_aggregation not in ("scalar_rbs", "region_hidden"):
        raise ValueError(f"unsupported gene_aggregation: {gene_aggregation!r}")
    training_cfg = config.get("training", {})
    age_loss_weight = float(training_cfg.get("age_loss_weight", 1.0))
    tissue_loss_weight = float(training_cfg.get("tissue_loss_weight", 1.0))
    sex_loss_weight = float(training_cfg.get("sex_loss_weight", 1.0))
    patience_raw = training_cfg.get("early_stopping_patience")
    early_stopping_patience = None if patience_raw is None else int(patience_raw)
    early_stopping_min_delta = float(training_cfg.get("early_stopping_min_delta", 0.0))
    lr = float(training_cfg.get("learning_rate", 1e-2))
    train_batch_size_raw = training_cfg.get("batch_size", "auto")
    train_batch_size: int | str | None
    if train_batch_size_raw is None:
        train_batch_size = "auto"
    elif isinstance(train_batch_size_raw, str):
        train_batch_size = train_batch_size_raw
    else:
        train_batch_size = int(train_batch_size_raw)
    gpu_share = max(1, int(training_cfg.get("gpu_share", 1)))
    fusion_cfg = config.get("fusion")
    if fusion_cfg is not None and not isinstance(fusion_cfg, dict):
        raise ValueError("config fusion must be a mapping")
    gene_linked_only = bool(training_cfg.get("gene_linked_only", False))
    gene_allocation_raw = model_cfg.get(
        "gene_allocation",
        training_cfg.get("gene_allocation", "legacy_nearest"),
    )
    gene_allocation = str(gene_allocation_raw)
    if gene_allocation not in ("explicit_only", "bounded_nearest", "legacy_nearest"):
        raise ValueError(f"unsupported gene_allocation: {gene_allocation!r}")
    max_nearest_raw = model_cfg.get(
        "max_nearest_gene_bp",
        training_cfg.get("max_nearest_gene_bp"),
    )
    max_nearest_gene_bp = None if max_nearest_raw is None else int(max_nearest_raw)
    primary_evaluation = str(training_cfg.get("primary_evaluation", "late_fusion"))
    if primary_evaluation not in ("late_fusion", "mbs_e2e"):
        raise ValueError(f"unsupported training.primary_evaluation: {primary_evaluation!r}")
    checkpoint_selection_mode = str(
        training_cfg.get("checkpoint_selection", TISSUE_PRIMARY_SELECTION)
    )
    if checkpoint_selection_mode not in (AGE_PRIMARY_SELECTION, TISSUE_PRIMARY_SELECTION):
        raise ValueError(
            f"unsupported training.checkpoint_selection: {checkpoint_selection_mode!r}"
        )
    extra_fusion_raw = training_cfg.get("extra_fusion_modes") or []
    extra_fusion_modes = tuple(str(m) for m in extra_fusion_raw)
    milestone = str(config.get("milestone", config.get("experiment", {}).get("name", "7F")))
    if "7g" in milestone.lower() or "7g" in run_id.lower():
        milestone_tag = "7G"
        default_report = project_root / "reports" / "inspection" / "stage0_7g_methylation_eval"
    else:
        milestone_tag = "7F"
        default_report = project_root / "reports" / "inspection" / "stage0_7f_rbs_gene_direct"
    report_rel = config.get("report_dir")
    if report_dir is None and report_rel:
        report_dir = Path(str(report_rel))
        if not report_dir.is_absolute():
            report_dir = project_root / report_dir
    if report_dir is None:
        report_dir = default_report

    folds_path = artifact_root / "splits" / split_id / "folds.json"
    if not folds_path.is_file():
        raise FileNotFoundError(f"frozen folds missing: {folds_path}")
    pack = load_frozen_folds(folds_path)
    folds = list(pack["folds"])
    if max_folds is not None:
        folds = folds[: int(max_folds)]

    print(f"[cascade] loading matrix/graph max_loci={max_loci}", flush=True)
    matrix_paths = matrix_store_paths(data_root / "canonical" / "matrices" / matrix_id)
    sample_index = read_sample_index(matrix_paths.sample_index_path)
    locus_index = read_locus_index(matrix_paths.locus_index_path)
    lr_edges, regions = load_graph_tables(data_root / "canonical" / "graphs" / graph_id)
    genes_path = data_root / "canonical" / "graphs" / graph_id / "genes.parquet"
    genes = pd.read_parquet(genes_path) if genes_path.is_file() else pd.DataFrame()
    print(f"[cascade] building assignment (regions={len(regions)})", flush=True)
    assignment = build_cascade_assignment(
        locus_index=locus_index,
        locus_region_edges=lr_edges,
        regions=regions,
        genes=genes,
        max_loci=max_loci,
        gene_allocation=gene_allocation,  # type: ignore[arg-type]
        max_nearest_gene_bp=max_nearest_gene_bp,
    )
    print(
        f"[cascade] assignment genes={assignment.n_genes} regions={assignment.n_regions} "
        f"orphan_rbs={assignment.n_orphan_rbs} direct={assignment.n_direct} "
        f"gene_allocation={gene_allocation}",
        flush=True,
    )
    betas_z = open_betas_zarr(matrix_paths.betas_path)
    row_by_id = {
        str(sid): int(row)
        for sid, row in zip(
            sample_index["sample_id"].astype(str),
            sample_index["row_index"].astype(int),
            strict=True,
        )
    }
    pheno_rel = Path(
        str(
            config.get("sample_phenotype_table")
            or pilot.get(
                "sample_phenotype_table",
                "canonical/phenotypes/sample_phenotype_table_age_tissue_sex_full_v1.parquet",
            )
        )
    )
    pheno_path = pheno_rel if pheno_rel.is_absolute() else data_root / pheno_rel
    phenotypes, class_names = load_multitask_phenotypes(pheno_path)
    ph_by_id = {p.sample_id: p for p in phenotypes}

    run_root = run_dir(artifact_root, run_id)
    run_root.mkdir(parents=True, exist_ok=True)
    fold_summaries: list[dict[str, Any]] = []
    n_cols = assignment.n_study_loci

    # Dense prefix load once (~3.5 GB float32 for 13.5k × 65k).
    print(f"[cascade] loading betas[:, :{n_cols}] into RAM…", flush=True)
    betas_all = np.asarray(betas_z[:, :n_cols], dtype=np.float32)
    print(f"[cascade] betas shape={betas_all.shape} dtype={betas_all.dtype}", flush=True)

    def _load_betas(ids: list[str]) -> tuple[np.ndarray, list[str]]:
        kept = [sid for sid in ids if sid in row_by_id]
        if not kept:
            raise ValueError("no samples with matrix rows")
        rows = np.asarray([row_by_id[sid] for sid in kept], dtype=np.int64)
        return betas_all[rows], kept

    for fold_i, fold in enumerate(folds):
        train_ids = [s for s in fold["train_sample_ids"] if s in row_by_id and s in ph_by_id]
        external_test_ids = fold.get("external_test_sample_ids") or []
        validation_ids = fold.get("validation_sample_ids") or []
        if external_test_ids:
            # Real held-out test set; validation stays a separate slice for
            # checkpoint selection (never used to pick the reported metrics).
            test_ids = [s for s in external_test_ids if s in row_by_id and s in ph_by_id]
            val_ids = [s for s in validation_ids if s in row_by_id and s in ph_by_id]
        else:
            # No external test in this fold; fall back to validation-as-test
            # (legacy behavior) and skip validation-based checkpointing.
            test_ids = [s for s in validation_ids if s in row_by_id and s in ph_by_id]
            val_ids = []
        if max_train_samples is not None and len(train_ids) > int(max_train_samples):
            rng = np.random.default_rng(seed + fold_i)
            pick = rng.choice(len(train_ids), size=int(max_train_samples), replace=False)
            train_ids = [train_ids[i] for i in sorted(pick.tolist())]
        if max_train_samples is not None and len(test_ids) > max(16, int(max_train_samples) // 4):
            n_te = max(16, int(max_train_samples) // 4)
            test_ids = test_ids[:n_te]
        if max_train_samples is not None and len(val_ids) > max(16, int(max_train_samples) // 4):
            n_va = max(16, int(max_train_samples) // 4)
            val_ids = val_ids[:n_va]
        print(
            f"[cascade] fold {fold_i} loading train={len(train_ids)} val={len(val_ids)} "
            f"test={len(test_ids)} loci={n_cols} epochs={max_epochs}",
            flush=True,
        )
        betas_tr, train_ids = _load_betas(train_ids)
        if val_ids:
            betas_va, val_ids = _load_betas(val_ids)
        else:
            betas_va = np.zeros((0, n_cols), dtype=betas_all.dtype)
        betas_te, test_ids = _load_betas(test_ids)
        betas = np.concatenate([betas_tr, betas_va, betas_te], axis=0)
        sample_ids = train_ids + val_ids + test_ids
        train_idx = np.arange(len(train_ids), dtype=np.int64)
        val_idx = np.arange(len(train_ids), len(train_ids) + len(val_ids), dtype=np.int64)
        test_idx = np.arange(len(train_ids) + len(val_ids), len(sample_ids), dtype=np.int64)

        ages = np.asarray(
            [float(ph_by_id[s].age or 0.0) for s in sample_ids], dtype=np.float64
        )
        tissue = np.asarray(
            [
                int(ph_by_id[s].class_index) if ph_by_id[s].tissue_mask else 0
                for s in sample_ids
            ],
            dtype=np.int64,
        )
        sex = np.asarray(
            [
                int(ph_by_id[s].sex_class_index or 0) if ph_by_id[s].sex_mask else 0
                for s in sample_ids
            ],
            dtype=np.int64,
        )
        age_mask = np.asarray([bool(ph_by_id[s].age_mask) for s in sample_ids], dtype=bool)
        tissue_mask = np.asarray(
            [bool(ph_by_id[s].tissue_mask) for s in sample_ids], dtype=bool
        )
        sex_mask = np.asarray([bool(ph_by_id[s].sex_mask) for s in sample_ids], dtype=bool)
        studies = np.asarray(
            [str(ph_by_id[s].study_id or "NA") for s in sample_ids],
            dtype=object,
        )
        fold_out = run_root / f"fold_{fold_i}"
        print(
            f"[cascade] fold {fold_i} train={len(train_ids)} val={len(val_ids)} "
            f"test={len(test_ids)} loci={n_cols} epochs={max_epochs}",
            flush=True,
        )
        metrics = train_cascade_on_arrays(
            assignment=assignment,
            betas=betas,
            train_idx=train_idx,
            test_idx=test_idx,
            ages=ages,
            tissue=tissue,
            sex=sex,
            study_ids=studies,
            sample_ids=sample_ids,
            class_names=list(class_names) if class_names else ["A", "B"],
            out_dir=fold_out,
            max_epochs=max_epochs,
            seed=seed + fold_i,
            device_str=device_str,
            lr=lr,
            age_mask=age_mask,
            tissue_mask=tissue_mask,
            sex_mask=sex_mask,
            cpg_hidden_dim=cpg_hidden,
            region_hidden_dim=region_hidden,
            dropout=dropout,
            cpg_pool=cpg_pool,
            region_pool=region_pool,
            gene_aggregation=gene_aggregation,
            gene_allocation_policy=gene_allocation,
            checkpoint_selection_mode=checkpoint_selection_mode,
            skip_if_done=skip_if_done,
            val_idx=val_idx,
            age_loss_weight=age_loss_weight,
            tissue_loss_weight=tissue_loss_weight,
            sex_loss_weight=sex_loss_weight,
            early_stopping_patience=early_stopping_patience,
            early_stopping_min_delta=early_stopping_min_delta,
            fusion=fusion_cfg,
            gene_linked_only=gene_linked_only,
            primary_evaluation=cast(PrimaryEvaluation, primary_evaluation),
            extra_fusion_modes=cast(tuple[FusionBlockMode, ...], extra_fusion_modes),
            locus_ids=locus_index["locus_id"].astype(str).tolist()[:n_cols],
            eval_only=eval_only,
            train_batch_size=train_batch_size,
            gpu_share=gpu_share,
            include_mbs_enet=bool(training_cfg.get("stage_a_include_mbs_enet", False)),
        )
        metrics["fold_id"] = fold.get("fold_id", fold_i)
        fold_summaries.append(metrics)

    report_dir.mkdir(parents=True, exist_ok=True)
    summary = {
        "milestone": milestone_tag,
        "topology": "rbs_gene_direct",
        "arm": "N-cascade-l1",
        "tbs_arm": False,
        "split_id": split_id,
        "matrix_id": matrix_id,
        "graph_id": graph_id,
        "max_loci": max_loci,
        "max_epochs": max_epochs,
        "n_restarts": 1,
        "remaining_ceiling": {
            "n_loci_in_matrix": int(locus_index.shape[0]),
            "n_loci_used": max_loci,
            "note": (
                "Prefix of matrix columns (not a claim that later CpGs are useless). "
                "Full-matrix train and a 2nd restart are deferred."
            ),
        },
        "encoder": {
            "cpg_hidden_dim": cpg_hidden,
            "region_hidden_dim": region_hidden,
            "dropout": dropout,
            "activation": "gelu",
            "layer_norm": True,
            "cpg_pool": cpg_pool,
            "region_pool": region_pool,
        },
        "training": {
            "age_loss_weight": age_loss_weight,
            "tissue_loss_weight": tissue_loss_weight,
            "sex_loss_weight": sex_loss_weight,
            "early_stopping_patience": early_stopping_patience,
            "early_stopping_min_delta": early_stopping_min_delta,
        },
        "assignment": {
            "n_genes": assignment.n_genes,
            "n_regions": assignment.n_regions,
            "n_orphan_rbs": assignment.n_orphan_rbs,
            "n_direct": assignment.n_direct,
        },
        "folds": fold_summaries,
        "adr": "0009",
    }
    write_json(report_dir / "summary.json", summary)
    (report_dir / "analysis.md").write_text(
        f"# Stage 0 Milestone {milestone_tag} — Hub cascade\n\n"
        f"Split `{split_id}`; max_loci={max_loci}; max_epochs={max_epochs}; "
        f"encoder {cpg_hidden}/{region_hidden}. "
        "No TBS arm. See `summary.json`.\n",
        encoding="utf-8",
    )
    score_dir = Path(str(fold_summaries[0]["score_dir"])) if fold_summaries else run_root
    return CascadeTrainResult(metrics=summary, score_dir=score_dir, report_dir=report_dir)


def train_cascade(
    *,
    project_root: Path,
    data_root: Path,
    artifact_root: Path,
    config: dict[str, Any] | None = None,
    config_path: Path | None = None,
    run_id: str = "stage0-7f-cascade",
    device_str: str = "cpu",
    overfit_fixture: bool = False,
    max_folds: int | None = None,
    max_train_samples: int | None = None,
    report_dir: Path | None = None,
    skip_if_done: bool = True,
    eval_only: bool = False,
) -> CascadeTrainResult:
    """CLI entry: fixture or Hub cascade on frozen 7E folds."""
    if overfit_fixture:
        return run_cascade_fixture(
            project_root=project_root,
            artifact_root=artifact_root,
            run_id=run_id,
            device_str=device_str,
        )
    if config is None:
        if config_path is None:
            raise ValueError("config or config_path required for Hub cascade")
        config = load_experiment_config(config_path)
    return run_cascade_hub(
        project_root=project_root,
        data_root=data_root,
        artifact_root=artifact_root,
        config=config,
        run_id=run_id,
        device_str=device_str,
        max_folds=max_folds,
        max_train_samples=max_train_samples,
        report_dir=report_dir,
        skip_if_done=skip_if_done,
        eval_only=eval_only,
    )
