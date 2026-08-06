"""Masked multitask heads and loss on shared FlatDeepSet MBS (Milestones 5c/5d)."""

from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn.functional as F
from torch import Tensor, nn

from mbs.models import SeedMaskedLinearHead
from mbs.training.dataset import FlatBatch


class MultitaskHeads(nn.Module):
    """DeepRVAT-style phenotype modules on shared gene MBS.

    Always instantiates parallel linear trait modules; missing phenotypes are
    handled only in the masked loss (not by switching heads per sample).
    ``PhenotypeModules`` is the preferred alias.
    """

    def __init__(
        self,
        n_genes: int,
        n_tissue_classes: int,
        *,
        n_sex_classes: int = 2,
        seed_mask: Tensor | None = None,
        sex_enabled: bool = False,
    ) -> None:
        super().__init__()
        self.n_genes = n_genes
        self.n_tissue_classes = n_tissue_classes
        self.n_sex_classes = int(n_sex_classes)
        self.sex_enabled = bool(sex_enabled)
        self.age_head = nn.Linear(n_genes, 1)
        mask = (
            seed_mask
            if seed_mask is not None
            else torch.ones(n_tissue_classes, n_genes, dtype=torch.float32)
        )
        self.tissue_head = SeedMaskedLinearHead(n_genes, n_tissue_classes, mask)
        self.sex_head = nn.Linear(n_genes, self.n_sex_classes) if self.sex_enabled else None
        with torch.no_grad():
            self.age_head.weight.normal_(0.0, 0.05)
            self.age_head.bias.zero_()
            self.tissue_head.gene_weight.normal_(0.0, 0.05)
            if self.sex_head is not None:
                self.sex_head.weight.normal_(0.0, 0.05)
                self.sex_head.bias.zero_()

    def forward_age(self, mbs: Tensor, present: Tensor) -> Tensor:
        masked = mbs * present.to(dtype=mbs.dtype)
        return self.age_head(masked).squeeze(-1)

    def forward_tissue(self, mbs: Tensor, present: Tensor) -> Tensor:
        return self.tissue_head(mbs, present)

    def forward_sex(self, mbs: Tensor, present: Tensor) -> Tensor:
        if self.sex_head is None:
            raise RuntimeError("sex head is disabled")
        masked = mbs * present.to(dtype=mbs.dtype)
        return self.sex_head(masked)


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
    lambda_sex: float = 1.0,
    huber_delta: float = 1.0,
    age_loss: str = "huber",
    class_weights: Tensor | None = None,
) -> MultitaskLossResult:
    """Factored DeepRVAT-style loss: sum only observed sample×trait terms.

    All phenotype modules run on the shared MBS; masks gate which terms enter
    the scalar loss (and thus which module parameters receive gradients).
    """
    if mbs.ndim == 1:
        mbs_b = mbs.unsqueeze(0)
        present_b = present.unsqueeze(0)
    else:
        mbs_b = mbs
        present_b = present

    device = mbs_b.device
    batch_size = int(mbs_b.shape[0])
    age_mask = batch.age_mask.reshape(-1).to(device=device, dtype=torch.bool)
    tissue_mask = batch.tissue_mask.reshape(-1).to(device=device, dtype=torch.bool)
    if batch.sex_mask is None:
        sex_mask = torch.zeros(batch_size, dtype=torch.bool, device=device)
    else:
        sex_mask = batch.sex_mask.reshape(-1).to(device=device, dtype=torch.bool)
    if (
        age_mask.numel() != batch_size
        or tissue_mask.numel() != batch_size
        or sex_mask.numel() != batch_size
    ):
        raise ValueError(
            f"mask length mismatch: batch={batch_size} "
            f"age_mask={age_mask.numel()} tissue_mask={tissue_mask.numel()} "
            f"sex_mask={sex_mask.numel()}"
        )

    total = torch.zeros((), device=device, dtype=torch.float32)
    metrics: dict[str, float] = {
        "age_loss": 0.0,
        "tissue_loss": 0.0,
        "sex_loss": 0.0,
        "mae": 0.0,
        "tissue_correct": 0.0,
        "sex_correct": 0.0,
        "age_n": 0.0,
        "tissue_n": 0.0,
        "sex_n": 0.0,
    }

    age_n = int(age_mask.sum().item())
    tissue_n = int(tissue_mask.sum().item())
    sex_n = int(sex_mask.sum().item())

    if age_n > 0:
        if batch.age_target is None:
            raise RuntimeError("age_mask set but age_target is None")
        pred = heads.forward_age(mbs_b, present_b)
        target = batch.age_target.reshape(-1).to(device=device, dtype=pred.dtype)
        pred_on = pred[age_mask]
        target_on = target[age_mask]
        if age_loss == "mse":
            age_term = F.mse_loss(pred_on, target_on)
        else:
            age_term = F.huber_loss(pred_on, target_on, delta=huber_delta)
        total = total + float(lambda_age) * age_term
        metrics["age_loss"] = float(age_term.detach().item())
        metrics["mae"] = float((pred_on.detach() - target_on.detach()).abs().mean().item())
        metrics["age_n"] = float(age_n)

    if tissue_n > 0:
        logits = heads.forward_tissue(mbs_b, present_b)
        targets = batch.tissue_target.reshape(-1).to(device=device)
        logits_on = logits[tissue_mask]
        targets_on = targets[tissue_mask]
        tissue_term = F.cross_entropy(
            logits_on,
            targets_on,
            weight=class_weights.to(device) if class_weights is not None else None,
        )
        total = total + float(lambda_tissue) * tissue_term
        metrics["tissue_loss"] = float(tissue_term.detach().item())
        pred_cls = logits_on.argmax(dim=-1)
        metrics["tissue_correct"] = float((pred_cls == targets_on).sum().item())
        metrics["tissue_n"] = float(tissue_n)

    if sex_n > 0:
        if not heads.sex_enabled:
            raise RuntimeError("sex_mask set but MultitaskHeads.sex_enabled is False")
        logits = heads.forward_sex(mbs_b, present_b)
        if batch.sex_target is None:
            raise RuntimeError("sex_mask set but sex_target is None")
        targets = batch.sex_target.reshape(-1).to(device=device)
        logits_on = logits[sex_mask]
        targets_on = targets[sex_mask]
        sex_term = F.cross_entropy(logits_on, targets_on)
        total = total + float(lambda_sex) * sex_term
        metrics["sex_loss"] = float(sex_term.detach().item())
        pred_cls = logits_on.argmax(dim=-1)
        metrics["sex_correct"] = float((pred_cls == targets_on).sum().item())
        metrics["sex_n"] = float(sex_n)

    if age_n == 0 and tissue_n == 0 and sex_n == 0:
        total = total + (mbs_b.sum() * 0.0)

    metrics["loss"] = float(total.detach().item())
    return MultitaskLossResult(loss=total, metrics=metrics)


# Preferred DeepRVAT-aligned names (aliases; keep Multitask* for 5c callers).
PhenotypeModules = MultitaskHeads
masked_phenotype_loss = masked_multitask_loss
PhenotypeLossResult = MultitaskLossResult
