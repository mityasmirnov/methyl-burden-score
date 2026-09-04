"""Masked multitask heads and loss on shared FlatDeepSet MBS (Milestones 5c/5d)."""

from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn.functional as F
from torch import Tensor, nn

from mbs.models import SeedMaskedLinearHead, center_mask_scores
from mbs.training.dataset import FlatBatch
from mbs.training.hier_dataset import HierBatch

# ADR 0011: a provided seed mask that selects too few genes is a config error,
# not a valid tiny panel. Fail closed below this many genes per trait row.
MIN_SEED_GENES = 32


def _prepare_seed_mask(mask: Tensor, *, n_outputs: int, n_genes: int, trait: str) -> Tensor:
    """Validate a provided seed mask (ADR 0011) and coerce to (n_outputs, n_genes)."""
    m = mask.detach().to(torch.float32)
    if m.ndim == 1:
        m = m.unsqueeze(0).expand(n_outputs, n_genes).contiguous()
    if tuple(m.shape) != (n_outputs, n_genes):
        raise ValueError(
            f"{trait} seed mask must have shape {(n_outputs, n_genes)}, found {tuple(m.shape)}"
        )
    row_sums = m.sum(dim=1)
    min_selected = int(row_sums.min().item()) if row_sums.numel() else 0
    if min_selected < MIN_SEED_GENES:
        raise ValueError(
            f"{trait} seed mask selects {min_selected} genes (< {MIN_SEED_GENES}); "
            "fail closed per ADR 0011"
        )
    return m


def _init_trait_head(head: SeedMaskedLinearHead | nn.Linear) -> None:
    """Symmetry-breaking init for either a masked or dense trait head."""
    if isinstance(head, SeedMaskedLinearHead):
        head.gene_weight.normal_(0.0, 0.05)
    else:
        head.weight.normal_(0.0, 0.05)
        head.bias.zero_()


class MultitaskHeads(nn.Module):
    """DeepRVAT-style phenotype modules on shared gene MBS.

    Always instantiates parallel linear trait modules; missing phenotypes are
    handled only in the masked loss (not by switching heads per sample).
    ``PhenotypeModules`` is the preferred alias.

    Age, tissue, and sex each accept an optional gene **seed mask** (ADR 0011).
    When a mask is provided the head becomes a :class:`SeedMaskedLinearHead` so
    only selected genes contribute gradients; when it is ``None`` age/sex keep a
    dense ``nn.Linear`` (backward compatible) and tissue keeps its all-ones mask.
    """

    def __init__(
        self,
        n_genes: int,
        n_tissue_classes: int,
        *,
        n_sex_classes: int = 2,
        seed_mask: Tensor | None = None,
        tissue_seed_mask: Tensor | None = None,
        age_seed_mask: Tensor | None = None,
        sex_seed_mask: Tensor | None = None,
        sex_enabled: bool = False,
        n_disease_labels: int = 0,
        n_cancer_labels: int = 0,
        neutral_score: float = 0.5,
    ) -> None:
        super().__init__()
        self.n_genes = n_genes
        self.n_tissue_classes = n_tissue_classes
        self.n_sex_classes = int(n_sex_classes)
        self.sex_enabled = bool(sex_enabled)
        self.n_disease_labels = int(n_disease_labels)
        self.n_cancer_labels = int(n_cancer_labels)
        self.neutral_score = float(neutral_score)

        # Age: masked head only when a mask is supplied; else dense (compat).
        self.age_head: SeedMaskedLinearHead | nn.Linear
        if age_seed_mask is not None:
            age_mask = _prepare_seed_mask(age_seed_mask, n_outputs=1, n_genes=n_genes, trait="age")
            self.age_head = SeedMaskedLinearHead(
                n_genes, 1, age_mask, neutral_score=self.neutral_score
            )
        else:
            self.age_head = nn.Linear(n_genes, 1)

        # Tissue: always a masked head; ``tissue_seed_mask`` wins over the
        # ``seed_mask`` alias, and an all-ones mask is the dense default.
        tissue_mask_arg = tissue_seed_mask if tissue_seed_mask is not None else seed_mask
        if tissue_mask_arg is not None:
            mask = _prepare_seed_mask(
                tissue_mask_arg, n_outputs=n_tissue_classes, n_genes=n_genes, trait="tissue"
            )
        else:
            mask = torch.ones(n_tissue_classes, n_genes, dtype=torch.float32)
        self.tissue_head = SeedMaskedLinearHead(
            n_genes, n_tissue_classes, mask, neutral_score=self.neutral_score
        )

        # Sex: masked head only when a mask is supplied; else dense (compat).
        self.sex_head: SeedMaskedLinearHead | nn.Linear | None
        if not self.sex_enabled:
            self.sex_head = None
        elif sex_seed_mask is not None:
            sex_mask = _prepare_seed_mask(
                sex_seed_mask, n_outputs=self.n_sex_classes, n_genes=n_genes, trait="sex"
            )
            self.sex_head = SeedMaskedLinearHead(
                n_genes, self.n_sex_classes, sex_mask, neutral_score=self.neutral_score
            )
        else:
            self.sex_head = nn.Linear(n_genes, self.n_sex_classes)

        self.disease_head = (
            nn.Linear(n_genes, self.n_disease_labels) if self.n_disease_labels > 0 else None
        )
        self.cancer_head = (
            nn.Linear(n_genes, self.n_cancer_labels) if self.n_cancer_labels > 0 else None
        )
        with torch.no_grad():
            _init_trait_head(self.age_head)
            self.tissue_head.gene_weight.normal_(0.0, 0.05)
            if self.sex_head is not None:
                _init_trait_head(self.sex_head)
            if self.disease_head is not None:
                self.disease_head.weight.normal_(0.0, 0.05)
                self.disease_head.bias.zero_()
            if self.cancer_head is not None:
                self.cancer_head.weight.normal_(0.0, 0.05)
                self.cancer_head.bias.zero_()

    def _centered(self, mbs: Tensor, present: Tensor) -> Tensor:
        return center_mask_scores(mbs, present, neutral_score=self.neutral_score)

    def forward_age(self, mbs: Tensor, present: Tensor) -> Tensor:
        if isinstance(self.age_head, SeedMaskedLinearHead):
            return self.age_head(mbs, present).squeeze(-1)
        return self.age_head(self._centered(mbs, present)).squeeze(-1)

    def forward_tissue(self, mbs: Tensor, present: Tensor) -> Tensor:
        return self.tissue_head(mbs, present)

    def forward_sex(self, mbs: Tensor, present: Tensor) -> Tensor:
        if self.sex_head is None:
            raise RuntimeError("sex head is disabled")
        if isinstance(self.sex_head, SeedMaskedLinearHead):
            return self.sex_head(mbs, present)
        return self.sex_head(self._centered(mbs, present))

    def forward_disease(self, mbs: Tensor, present: Tensor) -> Tensor:
        if self.disease_head is None:
            raise RuntimeError("disease head is disabled")
        return self.disease_head(self._centered(mbs, present))

    def forward_cancer(self, mbs: Tensor, present: Tensor) -> Tensor:
        if self.cancer_head is None:
            raise RuntimeError("cancer head is disabled")
        return self.cancer_head(self._centered(mbs, present))


@dataclass(frozen=True, slots=True)
class MultitaskLossResult:
    loss: Tensor
    metrics: dict[str, float]


def masked_multitask_loss(
    *,
    mbs: Tensor,
    present: Tensor,
    heads: MultitaskHeads,
    batch: FlatBatch | HierBatch,
    lambda_age: float = 1.0,
    lambda_tissue: float = 1.0,
    lambda_sex: float = 1.0,
    lambda_disease: float = 1.0,
    lambda_cancer: float = 1.0,
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
        "disease_n": 0.0,
        "cancer_n": 0.0,
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

    disease_n = 0
    cancer_n = 0
    disease_mask = batch.disease_mask
    if disease_mask is not None and heads.disease_head is not None:
        dmask = disease_mask.to(device=device, dtype=torch.bool)
        if dmask.any():
            if batch.disease_target is None:
                raise RuntimeError("disease_mask set but disease_target is None")
            logits = heads.forward_disease(mbs_b, present_b)
            target = batch.disease_target.to(device=device, dtype=logits.dtype)
            disease_term = F.binary_cross_entropy_with_logits(logits[dmask], target[dmask])
            total = total + float(lambda_disease) * disease_term
            metrics["disease_loss"] = float(disease_term.detach().item())
            disease_n = int(dmask.sum().item())
            metrics["disease_n"] = float(disease_n)
    cancer_mask = batch.cancer_mask
    if cancer_mask is not None and heads.cancer_head is not None:
        cmask = cancer_mask.to(device=device, dtype=torch.bool)
        if cmask.any():
            if batch.cancer_target is None:
                raise RuntimeError("cancer_mask set but cancer_target is None")
            logits = heads.forward_cancer(mbs_b, present_b)
            target = batch.cancer_target.to(device=device, dtype=logits.dtype)
            cancer_term = F.binary_cross_entropy_with_logits(logits[cmask], target[cmask])
            total = total + float(lambda_cancer) * cancer_term
            metrics["cancer_loss"] = float(cancer_term.detach().item())
            cancer_n = int(cmask.sum().item())
            metrics["cancer_n"] = float(cancer_n)

    if age_n == 0 and tissue_n == 0 and sex_n == 0 and disease_n == 0 and cancer_n == 0:
        total = total + (mbs_b.sum() * 0.0)

    metrics["loss"] = float(total.detach().item())
    return MultitaskLossResult(loss=total, metrics=metrics)


# Preferred DeepRVAT-aligned names (aliases; keep Multitask* for 5c callers).
PhenotypeModules = MultitaskHeads
masked_phenotype_loss = masked_multitask_loss
PhenotypeLossResult = MultitaskLossResult
