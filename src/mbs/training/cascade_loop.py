"""Milestone 7F cascade trainer: RBS→gene MBS + orphan RBS + leftover direct."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

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
from mbs.training.cascade_assign import CascadeAssignment, build_cascade_assignment
from mbs.training.cascade_scores import (
    fusion_feature_matrix,
    load_cascade_score_blocks,
    write_cascade_score_dir,
)
from mbs.training.dev_cv import DEFAULT_SPLIT_ID, load_frozen_folds
from mbs.training.direct_cpg import direct_cpg_design_matrix, fit_direct_elasticnet
from mbs.training.features import beta_to_m_value
from mbs.training.late_fusion import evaluate_late_fusion
from mbs.training.locus_gene import load_graph_tables
from mbs.training.loop import load_experiment_config, resolve_device
from mbs.training.multitask import MultitaskHeads
from mbs.training.phenotypes import load_multitask_phenotypes
from mbs.training.run_artifacts import run_dir


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


def _forward_sample(
    model: CascadeDeepSet,
    assignment: CascadeAssignment,
    beta_row: np.ndarray,
    *,
    device: torch.device,
) -> dict[str, torch.Tensor]:
    feats = _dense_cpg_features(beta_row)
    cols = assignment.edge_col_index
    if cols.size == 0:
        cpg_features = torch.zeros(0, feats.shape[1], device=device, dtype=torch.float32)
        cpg_to_region = torch.zeros(0, device=device, dtype=torch.long)
    else:
        cpg_features = torch.from_numpy(feats[cols]).to(device)
        cpg_to_region = torch.from_numpy(assignment.edge_region_index.astype(np.int64)).to(device)
    region_type = torch.from_numpy(assignment.region_type_id.astype(np.int64)).to(device)
    region_to_gene = torch.from_numpy(assignment.region_to_gene.astype(np.int64)).to(device)
    orphan_idx = torch.from_numpy(assignment.orphan_region_indices).to(device)
    return model(
        cpg_features=cpg_features,
        cpg_to_region=cpg_to_region,
        region_type=region_type,
        region_to_gene=region_to_gene,
        n_regions=assignment.n_regions,
        n_gene_instances=max(assignment.n_genes, 1),
        orphan_region_indices=orphan_idx,
    )


def score_samples(
    model: CascadeDeepSet,
    assignment: CascadeAssignment,
    betas: np.ndarray,
    *,
    device: torch.device,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return mbs, gene_present, orphan_rbs for all sample rows."""
    model.eval()
    n = betas.shape[0]
    n_genes = max(assignment.n_genes, 1)
    n_orphan = assignment.n_orphan_rbs
    mbs_out = np.full((n, n_genes), 0.5, dtype=np.float32)
    present_out = np.zeros((n, n_genes), dtype=bool)
    orphan_out = np.full((n, n_orphan), 0.5, dtype=np.float32)
    with torch.no_grad():
        for i in range(n):
            out = _forward_sample(model, assignment, betas[i], device=device)
            m = out["mbs"].detach().cpu().numpy().astype(np.float32)
            p = out["present"].detach().cpu().numpy().astype(bool)
            if assignment.n_genes == 0:
                mbs_out[i, 0] = 0.5
                present_out[i, 0] = False
            else:
                mbs_out[i, : assignment.n_genes] = m[: assignment.n_genes]
                present_out[i, : assignment.n_genes] = p[: assignment.n_genes]
            if n_orphan:
                orphan_out[i] = out["orphan_rbs"].detach().cpu().numpy().astype(np.float32)
    if assignment.n_genes == 0:
        # Keep a single neutral gene column for head plumbing.
        return mbs_out, present_out, orphan_out
    return mbs_out[:, : assignment.n_genes], present_out[:, : assignment.n_genes], orphan_out


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
) -> dict[str, Any]:
    """Cheap proxy metrics from the model's own heads on a held-out validation slice."""
    from mbs.evaluation.metrics import multiclass_metrics, regression_metrics  # noqa: PLC0415

    out: dict[str, Any] = {"tissue_macro_f1": None, "age_mae": None}
    if betas_val.shape[0] == 0:
        return out
    mbs_v, present_v, _ = score_samples(model, assignment, betas_val, device=device)
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
    return out


def _validation_rank(val_metrics: dict[str, Any]) -> tuple[float, float]:
    """Higher is better: (tissue macro-F1, -age MAE), missing -> worst."""
    f1 = val_metrics.get("tissue_macro_f1")
    mae = val_metrics.get("age_mae")
    return (
        float(f1) if f1 is not None else -1.0,
        -float(mae) if mae is not None else -1e9,
    )


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
    cpg_hidden_dim: int = 64,
    region_hidden_dim: int = 32,
    dropout: float = 0.1,
    skip_if_done: bool = False,
    val_idx: np.ndarray | None = None,
    age_loss_weight: float = 1.0,
    tissue_loss_weight: float = 1.0,
    sex_loss_weight: float = 1.0,
    fusion: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Train CascadeDeepSet + MBS heads; write scores; late-fuse; return metrics."""
    score_dir = out_dir / "scores"
    manifest_path = score_dir / "score_manifest.json"
    metrics_path = out_dir / "metrics.json"
    if skip_if_done and manifest_path.is_file() and metrics_path.is_file():
        cached = json.loads(metrics_path.read_text(encoding="utf-8"))
        cached["skipped"] = True
        cached["score_dir"] = str(score_dir)
        return cached

    _set_seed(seed)
    device = resolve_device(device_str, require_cuda=False)
    n_region_types = max(len(assignment.region_types), 1)
    n_genes = max(assignment.n_genes, 1)
    model = CascadeDeepSet(
        1,
        n_region_types,
        cpg_hidden_dim=int(cpg_hidden_dim),
        region_hidden_dim=int(region_hidden_dim),
        dropout=float(dropout),
        activation="gelu",
        layer_norm=True,
    )
    heads = MultitaskHeads(n_genes, max(len(class_names), 2), sex_enabled=True)
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

    val_idx_a = None if val_idx is None else np.asarray(val_idx, dtype=np.int64)
    has_val = val_idx_a is not None and val_idx_a.size > 0
    out_dir.mkdir(parents=True, exist_ok=True)
    ckpt_path = out_dir / "best.pt"

    def _save_checkpoint() -> None:
        torch.save(
            {
                "model": model.state_dict(),
                "heads": heads.state_dict(),
                "seed": seed,
                "max_epochs": max_epochs,
                "cpg_hidden_dim": cpg_hidden_dim,
                "region_hidden_dim": region_hidden_dim,
            },
            ckpt_path,
        )

    best_rank: tuple[float, float] | None = None
    best_epoch = -1
    val_history: list[dict[str, Any]] = []

    for _epoch in range(max_epochs):
        model.train()
        heads.train()
        order = train_idx.copy()
        np.random.shuffle(order)
        for i in order.tolist():
            if not (bool(age_mask_a[i]) or bool(tissue_mask_a[i]) or bool(sex_mask_a[i])):
                continue
            out = _forward_sample(model, assignment, betas[i], device=device)
            mbs = out["mbs"].unsqueeze(0)
            present = out["present"].unsqueeze(0)
            if assignment.n_genes == 0:
                mbs = torch.full((1, 1), 0.5, device=device)
                present = torch.zeros(1, 1, dtype=torch.bool, device=device)
            loss = torch.zeros((), device=device)
            if bool(age_mask_a[i]):
                age_t = torch.tensor([float(ages[i])], device=device)
                age_pred = heads.forward_age(mbs, present)
                loss = loss + age_loss_weight * F.huber_loss(age_pred, age_t)
            if bool(tissue_mask_a[i]):
                tissue_t = torch.tensor([int(tissue[i])], device=device, dtype=torch.long)
                tissue_pred = heads.forward_tissue(mbs, present)
                loss = loss + tissue_loss_weight * F.cross_entropy(
                    tissue_pred, tissue_t, weight=tissue_class_weights
                )
            if bool(sex_mask_a[i]):
                sex_t = torch.tensor([int(sex[i])], device=device, dtype=torch.long)
                sex_pred = heads.forward_sex(mbs, present)
                if sex_pred is not None:
                    loss = loss + sex_loss_weight * F.cross_entropy(sex_pred, sex_t)
            opt.zero_grad(set_to_none=True)
            loss.backward()
            opt.step()
        if (_epoch + 1) % max(1, max_epochs // 5) == 0 or _epoch == 0:
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
            )
            rank = _validation_rank(val_metrics)
            val_history.append({"epoch": _epoch + 1, "rank": list(rank), **val_metrics})
            if best_rank is None or rank > best_rank:
                best_rank = rank
                best_epoch = _epoch + 1
                _save_checkpoint()

    checkpoint_selection: dict[str, Any] = {
        "has_validation": has_val,
        "n_val": int(val_idx_a.size) if val_idx_a is not None else 0,
        "max_epochs": max_epochs,
        "val_history": val_history,
    }
    if has_val and best_epoch > 0:
        checkpoint_selection["best_epoch"] = best_epoch
        checkpoint_selection["selection"] = "validation_tissue_macro_f1_then_age_mae"
        # Reload the best-validation checkpoint (may not be the final epoch).
        ckpt = torch.load(ckpt_path, map_location=device)
        model.load_state_dict(ckpt["model"])
        heads.load_state_dict(ckpt["heads"])
    else:
        checkpoint_selection["best_epoch"] = max_epochs
        checkpoint_selection["selection"] = "final_epoch_no_validation"
        _save_checkpoint()

    mbs_all, present_all, orphan_all = score_samples(model, assignment, betas, device=device)
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
    )

    blocks = load_cascade_score_blocks(score_dir)
    x = fusion_feature_matrix(blocks)
    if "tbs" in blocks:
        raise RuntimeError("TBS leaked into 7F fusion")
    x_tr = x[train_idx]
    x_te = x[test_idx]
    fused = evaluate_late_fusion(
        scores_train=x_tr,
        scores_test=x_te,
        age_train=ages[train_idx],
        age_mask_train=age_mask_a[train_idx],
        tissue_train=tissue[train_idx],
        tissue_mask_train=tissue_mask_a[train_idx],
        sex_train=sex[train_idx],
        sex_mask_train=sex_mask_a[train_idx],
        age_test=ages[test_idx],
        age_mask_test=age_mask_a[test_idx],
        tissue_test=tissue[test_idx],
        tissue_mask_test=tissue_mask_a[test_idx],
        sex_test=sex[test_idx],
        sex_mask_test=sex_mask_a[test_idx],
        study_ids_test=study_ids[test_idx],
        tissue_class_names=list(class_names) if class_names else None,
        fusion=fusion,
    )
    fused["n_orphan_rbs"] = int(assignment.n_orphan_rbs)
    fused["n_direct"] = int(assignment.n_direct)
    fused["n_genes"] = int(assignment.n_genes)
    fused["tbs_arm"] = False
    fused["score_dir"] = str(score_dir)
    fused["fusion_n_features"] = int(x.shape[1])
    fused["skipped"] = False
    fused["checkpoint"] = str(ckpt_path)
    fused["checkpoint_selection"] = checkpoint_selection
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
    training_cfg = config.get("training", {})
    age_loss_weight = float(training_cfg.get("age_loss_weight", 1.0))
    tissue_loss_weight = float(training_cfg.get("tissue_loss_weight", 1.0))
    sex_loss_weight = float(training_cfg.get("sex_loss_weight", 1.0))
    lr = float(training_cfg.get("learning_rate", 1e-2))
    fusion_cfg = config.get("fusion")
    if fusion_cfg is not None and not isinstance(fusion_cfg, dict):
        raise ValueError("config fusion must be a mapping")
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
    )
    print(
        f"[cascade] assignment genes={assignment.n_genes} regions={assignment.n_regions} "
        f"orphan_rbs={assignment.n_orphan_rbs} direct={assignment.n_direct}",
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
            skip_if_done=skip_if_done,
            val_idx=val_idx,
            age_loss_weight=age_loss_weight,
            tissue_loss_weight=tissue_loss_weight,
            sex_loss_weight=sex_loss_weight,
            fusion=fusion_cfg,
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
    )
