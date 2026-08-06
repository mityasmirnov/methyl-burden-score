"""Unit tests for masked phenotype-module loss (Milestones 5c/5d)."""

from __future__ import annotations

import torch

from mbs.training.dataset import FlatBatch
from mbs.training.multitask import MultitaskHeads, masked_multitask_loss


def _batch(
    *,
    age_on: bool,
    tissue_on: bool,
    age: float = 0.0,
    tissue_cls: int = 0,
    n_genes: int = 4,
    sex_on: bool = False,
    sex_cls: int = 0,
) -> FlatBatch:
    return FlatBatch(
        sample_ids=["s0"],
        cpg_features=torch.zeros(1, 2),
        cpg_to_gene=torch.zeros(1, dtype=torch.long),
        n_genes=n_genes,
        tissue_target=torch.tensor([tissue_cls], dtype=torch.long),
        tissue_mask=torch.tensor([tissue_on]),
        age_target=torch.tensor([age], dtype=torch.float32) if age_on else None,
        age_mask=torch.tensor([age_on]),
        sex_target=torch.tensor([sex_cls], dtype=torch.long),
        sex_mask=torch.tensor([sex_on]),
    )


def test_masked_loss_ignores_unlabeled_age_head() -> None:
    n_genes = 4
    n_classes = 3
    heads = MultitaskHeads(n_genes, n_classes)
    mbs = torch.full((n_genes,), 0.7)
    present = torch.ones(n_genes, dtype=torch.bool)
    batch = _batch(age_on=False, tissue_on=True, tissue_cls=1, n_genes=n_genes)
    result = masked_multitask_loss(
        mbs=mbs,
        present=present,
        heads=heads,
        batch=batch,
        lambda_age=1.0,
        lambda_tissue=1.0,
    )
    assert result.metrics["age_n"] == 0.0
    assert result.metrics["tissue_n"] == 1.0
    assert result.metrics["age_loss"] == 0.0
    assert result.metrics["tissue_loss"] > 0.0
    assert torch.isfinite(result.loss)


def test_masked_loss_both_heads_when_masked() -> None:
    n_genes = 4
    heads = MultitaskHeads(n_genes, 3)
    mbs = torch.linspace(0.2, 0.8, n_genes)
    present = torch.ones(n_genes, dtype=torch.bool)
    batch = _batch(age_on=True, tissue_on=True, age=0.5, tissue_cls=0, n_genes=n_genes)
    result = masked_multitask_loss(
        mbs=mbs,
        present=present,
        heads=heads,
        batch=batch,
        lambda_age=1.0,
        lambda_tissue=1.0,
    )
    assert result.metrics["age_n"] == 1.0
    assert result.metrics["tissue_n"] == 1.0
    assert result.metrics["age_loss"] > 0.0
    assert result.metrics["tissue_loss"] > 0.0
    # Total should be sum of terms (lambdas=1).
    expected = result.metrics["age_loss"] + result.metrics["tissue_loss"]
    assert abs(result.metrics["loss"] - expected) < 1e-5


def test_zero_mask_batch_does_not_nan() -> None:
    n_genes = 4
    heads = MultitaskHeads(n_genes, 2)
    mbs = torch.randn(n_genes, requires_grad=True)
    present = torch.ones(n_genes, dtype=torch.bool)
    batch = _batch(age_on=False, tissue_on=False, n_genes=n_genes)
    result = masked_multitask_loss(
        mbs=mbs,
        present=present,
        heads=heads,
        batch=batch,
    )
    assert torch.isfinite(result.loss)
    result.loss.backward()
    assert mbs.grad is not None


def test_masked_loss_batched_partial_masks() -> None:
    n_genes = 4
    heads = MultitaskHeads(n_genes, 3)
    mbs = torch.randn(3, n_genes)
    present = torch.ones(3, n_genes, dtype=torch.bool)
    batch = FlatBatch(
        sample_ids=["a", "b", "c"],
        cpg_features=torch.zeros(3, 2),
        cpg_to_gene=torch.arange(3, dtype=torch.long),
        n_genes=n_genes,
        tissue_target=torch.tensor([0, 1, 2], dtype=torch.long),
        tissue_mask=torch.tensor([True, False, True]),
        age_target=torch.tensor([0.1, 0.2, 0.3], dtype=torch.float32),
        age_mask=torch.tensor([True, True, False]),
    )
    result = masked_multitask_loss(
        mbs=mbs,
        present=present,
        heads=heads,
        batch=batch,
    )
    assert result.metrics["age_n"] == 2.0
    assert result.metrics["tissue_n"] == 2.0
    assert torch.isfinite(result.loss)


def test_masked_sex_loss_only_when_observed() -> None:
    n_genes = 4
    heads = MultitaskHeads(n_genes, 3, sex_enabled=True)
    mbs = torch.full((n_genes,), 0.6)
    present = torch.ones(n_genes, dtype=torch.bool)
    off = _batch(age_on=False, tissue_on=False, sex_on=False, n_genes=n_genes)
    on = _batch(age_on=False, tissue_on=False, sex_on=True, sex_cls=1, n_genes=n_genes)
    r_off = masked_multitask_loss(mbs=mbs, present=present, heads=heads, batch=off)
    r_on = masked_multitask_loss(mbs=mbs, present=present, heads=heads, batch=on)
    assert r_off.metrics["sex_n"] == 0.0
    assert r_off.metrics["sex_loss"] == 0.0
    assert r_on.metrics["sex_n"] == 1.0
    assert r_on.metrics["sex_loss"] > 0.0
    assert torch.isfinite(r_on.loss)


def test_sex_module_no_grad_when_masked() -> None:
    n_genes = 4
    heads = MultitaskHeads(n_genes, 2, sex_enabled=True)
    mbs = torch.randn(n_genes, requires_grad=True)
    present = torch.ones(n_genes, dtype=torch.bool)
    batch = _batch(age_on=True, tissue_on=False, sex_on=False, age=0.0, n_genes=n_genes)
    result = masked_multitask_loss(mbs=mbs, present=present, heads=heads, batch=batch)
    result.loss.backward()
    assert heads.sex_head is not None
    assert heads.sex_head.weight.grad is None or torch.all(heads.sex_head.weight.grad == 0)
