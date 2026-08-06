"""Masked multitask heads and loss on shared FlatDeepSet MBS (Milestone 5c)."""

from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn.functional as F
from torch import Tensor, nn

from mbs.models import SeedMaskedLinearHead
from mbs.training.dataset import FlatBatch


class MultitaskHeads(nn.Module):
    """Linear age head + SeedMaskedLinearHead tissue CE on shared gene MBS."""

    def __init__(
        self,
        n_genes: int,
        n_tissue_classes: int,
        *,
        seed_mask: Tensor | None = None,
    ) -> None:
        super().__init__()
        self.n_genes = n_genes
        self.n_tissue_classes = n_tissue_classes
        self.age_head = nn.Linear(n_genes, 1)
        mask = (
            seed_mask
            if seed_mask is not None
            else torch.ones(n_tissue_classes, n_genes, dtype=torch.float32)
        )
        self.tissue_head = SeedMaskedLinearHead(n_genes, n_tissue_classes, mask)
        with torch.no_grad():
            self.age_head.weight.normal_(0.0, 0.05)
            self.age_head.bias.zero_()
            self.tissue_head.gene_weight.normal_(0.0, 0.05)

    def forward_age(self, mbs: Tensor, present: Tensor) -> Tensor:
        masked = mbs * present
        return self.age_head(masked).squeeze(-1)

    def forward_tissue(self, mbs: Tensor, present: Tensor) -> Tensor:
        return self.tissue_head(mbs, present)


@dataclass(frozen=True, slots=True)
class MultitaskLossResult:
    loss: Tensor
    metrics: dict[str, float]


def masked_multitask_loss(
    *,
    mbs: Tensor,
    present: Tensor,
    heads: MultitaskHeads,
    batch: FlatBatch,
    lambda_age: float = 1.0,
    lambda_tissue: float = 1.0,
    huber_delta: float = 1.0,
    age_loss: str = "huber",
    class_weights: Tensor | None = None,
) -> MultitaskLossResult:
    """Weighted sum of age + tissue losses; unlabeled heads contribute 0.

    ``mbs`` / ``present`` are shape ``[n_genes]`` (single-sample) or
    ``[batch, n_genes]``. Batch size 1 is the Stage 0 flat loop default.
    """
    if mbs.ndim == 1:
        mbs_b = mbs.unsqueeze(0)
        present_b = present.unsqueeze(0)
    else:
        mbs_b = mbs
        present_b = present

    device = mbs_b.device
    total = torch.zeros((), device=device, dtype=mbs_b.dtype)
    metrics: dict[str, float] = {
        "age_loss": 0.0,
        "tissue_loss": 0.0,
        "mae": 0.0,
        "tissue_correct": 0.0,
        "age_n": 0.0,
        "tissue_n": 0.0,
    }

    age_on = bool(batch.age_mask.reshape(-1)[0].item()) and batch.age_target is not None
    tissue_on = bool(batch.tissue_mask.reshape(-1)[0].item())

    if age_on:
        age_target = batch.age_target
        if age_target is None:
            raise RuntimeError("age_mask set but age_target is None")
        pred = heads.forward_age(mbs_b, present_b)
        target = age_target.reshape(pred.shape)
        if age_loss == "mse":
            age_term = F.mse_loss(pred, target)
        else:
            age_term = F.huber_loss(pred, target, delta=huber_delta)
        total = total + float(lambda_age) * age_term
        metrics["age_loss"] = float(age_term.detach().item())
        metrics["mae"] = float((pred.detach() - target.detach()).abs().mean().item())
        metrics["age_n"] = 1.0

    if tissue_on:
        logits = heads.forward_tissue(mbs_b, present_b)
        tissue_term = F.cross_entropy(
            logits,
            batch.tissue_target,
            weight=class_weights.to(device) if class_weights is not None else None,
        )
        total = total + float(lambda_tissue) * tissue_term
        metrics["tissue_loss"] = float(tissue_term.detach().item())
        pred_cls = int(logits.argmax(dim=-1).item())
        metrics["tissue_correct"] = float(pred_cls == int(batch.tissue_target.item()))
        metrics["tissue_n"] = 1.0

    if not age_on and not tissue_on:
        # Keep graph connected without contributing a real objective.
        total = total + (mbs_b.sum() * 0.0)

    metrics["loss"] = float(total.detach().item())
    return MultitaskLossResult(loss=total, metrics=metrics)
