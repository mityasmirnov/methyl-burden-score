"""7G′ Stage B flat-region fold trainer (FlatDeepSetRegion / N-light-type)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import torch
import torch.nn.functional as F

from mbs.training.transparent_baselines import evaluate_multitask_predictions
from mbs.models import FlatDeepSetRegion
from mbs.training.cascade_assign import CascadeAssignment, assignment_col_subset
from mbs.training.cascade_loop import _tissue_class_weights
from mbs.training.flat_region_features import (
    build_flat_region_gene_index,
    gather_flat_region_features,
)
from mbs.training.loop import resolve_device
from mbs.training.multitask import MultitaskHeads


def _set_seed(seed: int) -> None:
    torch.manual_seed(seed)
    np.random.seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def train_flat_region_on_arrays(
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
    max_epochs: int = 15,
    seed: int = 42,
    device_str: str = "cpu",
    lr: float = 1e-3,
    age_mask: np.ndarray | None = None,
    tissue_mask: np.ndarray | None = None,
    sex_mask: np.ndarray | None = None,
    panel_cols: np.ndarray | None = None,
    tissue_loss_weight: float = 3.0,
    age_loss_weight: float = 0.3,
    sex_loss_weight: float = 1.0,
) -> dict[str, Any]:
    """Train FlatDeepSetRegion on one fold; return test metrics."""
    out_dir.mkdir(parents=True, exist_ok=True)
    if panel_cols is not None:
        assignment = assignment_col_subset(assignment, panel_cols)
    index = build_flat_region_gene_index(assignment)
    n_genes = max(index.n_genes, 1)
    n_classes = max(len(class_names), 2)
    device = resolve_device(device_str, require_cuda=False)
    input_dim = 1 + len(index.region_types) + 1
    model = FlatDeepSetRegion(
        input_dim,
        phi_hidden_dim=64,
        phi_layers=2,
        rho_hidden_dim=10,
        rho_layers=2,
        pool="max",
    ).to(device)
    heads = MultitaskHeads(n_genes, n_classes, sex_enabled=True).to(device)
    opt = torch.optim.AdamW(list(model.parameters()) + list(heads.parameters()), lr=lr)
    age_mask_a = np.ones(len(sample_ids), dtype=bool) if age_mask is None else age_mask
    tissue_mask_a = np.ones(len(sample_ids), dtype=bool) if tissue_mask is None else tissue_mask
    sex_mask_a = np.zeros(len(sample_ids), dtype=bool) if sex_mask is None else sex_mask
    tissue_class_weights = _tissue_class_weights(
        tissue[train_idx],
        tissue_mask_a[train_idx],
        n_classes=n_classes,
        device=device,
    )

    def forward_sample(row: int) -> tuple[torch.Tensor, torch.Tensor]:
        feats, cpg_to_gene = gather_flat_region_features(
            beta_row=betas[row],
            index=index,
        )
        if feats.shape[0] == 0:
            mbs = torch.full((n_genes,), 0.5, device=device)
            present = torch.zeros(n_genes, dtype=torch.bool, device=device)
            return mbs.unsqueeze(0), present.unsqueeze(0)
        x = torch.from_numpy(feats).to(device)
        g = torch.from_numpy(cpg_to_gene).to(device)
        out = model(x, g, n_genes)
        return out["mbs"].unsqueeze(0), out["present"].unsqueeze(0)

    _set_seed(seed)
    train_idx = np.asarray(train_idx, dtype=np.int64)
    test_idx = np.asarray(test_idx, dtype=np.int64)
    for _epoch in range(max_epochs):
        model.train()
        heads.train()
        order = train_idx.copy()
        np.random.shuffle(order)
        for i in order.tolist():
            if not (bool(age_mask_a[i]) or bool(tissue_mask_a[i]) or bool(sex_mask_a[i])):
                continue
            mbs, present = forward_sample(int(i))
            loss = torch.zeros((), device=device)
            if bool(age_mask_a[i]):
                age_t = torch.tensor([float(ages[i])], device=device)
                loss = loss + age_loss_weight * F.huber_loss(heads.forward_age(mbs, present), age_t)
            if bool(tissue_mask_a[i]):
                tissue_t = torch.tensor([int(tissue[i])], device=device, dtype=torch.long)
                loss = loss + tissue_loss_weight * F.cross_entropy(
                    heads.forward_tissue(mbs, present),
                    tissue_t,
                    weight=tissue_class_weights,
                )
            if bool(sex_mask_a[i]):
                sex_t = torch.tensor([int(sex[i])], device=device, dtype=torch.long)
                sex_pred = heads.forward_sex(mbs, present)
                if sex_pred is not None:
                    loss = loss + sex_loss_weight * F.cross_entropy(sex_pred, sex_t)
            opt.zero_grad(set_to_none=True)
            loss.backward()
            opt.step()

    model.eval()
    heads.eval()
    tissue_pred: list[int] = []
    with torch.no_grad():
        for i in test_idx.tolist():
            mbs, present = forward_sample(int(i))
            logits = heads.forward_tissue(mbs, present)
            tissue_pred.append(int(logits.argmax(dim=-1).item()))
    test_idx_a = np.asarray(test_idx, dtype=np.int64)
    train_idx_a = np.asarray(train_idx, dtype=np.int64)
    tm_tr = tissue_mask_a[train_idx_a]
    tissue_valid_classes = (
        set(tissue[train_idx_a][tm_tr].tolist()) if tm_tr.any() else None
    )
    metrics = evaluate_multitask_predictions(
        preds={"tissue": np.asarray(tissue_pred, dtype=np.int64)},
        age=ages[test_idx_a],
        age_mask=age_mask_a[test_idx_a],
        tissue=tissue[test_idx_a],
        tissue_mask=tissue_mask_a[test_idx_a],
        sex=sex[test_idx_a],
        sex_mask=sex_mask_a[test_idx_a],
        study_ids=study_ids[test_idx_a],
        tissue_class_names=list(class_names),
        tissue_valid_classes=tissue_valid_classes,
    )
    payload = {
        "metrics": metrics,
        "arm": "N-light-type",
        "n_genes": index.n_genes,
        "eval_split": "test",
        "n_eval_samples": int(test_idx_a.size),
    }
    (out_dir / "metrics.json").write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    return payload
