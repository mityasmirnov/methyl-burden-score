"""Seed-mask behaviour for MultitaskHeads (ADR 0011)."""

from __future__ import annotations

import pytest
import torch

from mbs.training.checkpoint_selection import validation_rank, validation_rank_age_primary
from mbs.training.multitask import MIN_SEED_GENES, MultitaskHeads


def test_undersized_age_mask_raises() -> None:
    n_genes = 64
    mask = torch.zeros(1, n_genes)
    mask[:, :10] = 1.0
    with pytest.raises(ValueError, match="age seed mask"):
        MultitaskHeads(n_genes, 3, age_seed_mask=mask, sex_enabled=True)


def test_age_mask_zeros_unselected_gradient() -> None:
    n_genes = 64
    mask = torch.zeros(1, n_genes)
    mask[:, :MIN_SEED_GENES] = 1.0
    heads = MultitaskHeads(n_genes, 2, age_seed_mask=mask, sex_enabled=False)
    mbs = torch.randn(4, n_genes, requires_grad=True)
    present = torch.ones(4, n_genes)
    pred = heads.forward_age(mbs, present)
    pred.sum().backward()
    assert mbs.grad is not None
    # Unselected genes must receive ~0 gradient from the age head.
    assert float(mbs.grad[:, MIN_SEED_GENES:].abs().max()) < 1e-8
    assert float(mbs.grad[:, :MIN_SEED_GENES].abs().sum()) > 0.0


def test_validation_rank_age_primary_prefers_lower_mae() -> None:
    worse = {"mae": 12.0, "macro_f1": 0.40, "sex_auroc": 0.80}
    better = {"mae": 10.0, "macro_f1": 0.35, "sex_auroc": 0.70}
    assert validation_rank_age_primary(better) > validation_rank_age_primary(worse)
    # Legacy tissue-primary still prefers higher F1.
    assert validation_rank(worse) > validation_rank(better)
