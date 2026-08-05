from __future__ import annotations

import torch

from mbs.segment_ops import (
    segment_max,
    segment_mean,
    segment_softmax,
    segment_sqrt_sum,
    segment_sum,
)


def test_segment_reductions_and_empty_segments() -> None:
    values = torch.tensor(
        [
            [1.0, 2.0],
            [3.0, 0.0],
            [2.0, 5.0],
        ]
    )
    index = torch.tensor([0, 0, 2])

    total = segment_sum(values, index, n_segments=4)
    mean = segment_mean(values, index, n_segments=4)
    sqrt_sum = segment_sqrt_sum(values, index, n_segments=4)
    maximum, present = segment_max(values, index, n_segments=4)

    assert torch.allclose(
        total,
        torch.tensor(
            [
                [4.0, 2.0],
                [0.0, 0.0],
                [2.0, 5.0],
                [0.0, 0.0],
            ]
        ),
    )
    assert torch.allclose(mean[0], torch.tensor([2.0, 1.0]))
    assert torch.allclose(mean[1], torch.zeros(2))
    assert torch.allclose(sqrt_sum[0], torch.tensor([4.0, 2.0]) / torch.sqrt(torch.tensor(2.0)))
    assert torch.allclose(maximum[0], torch.tensor([3.0, 2.0]))
    assert torch.allclose(maximum[2], torch.tensor([2.0, 5.0]))
    assert present.tolist() == [True, False, True, False]


def test_permutation_invariance() -> None:
    generator = torch.Generator().manual_seed(7)
    values = torch.randn(20, 5, generator=generator)
    index = torch.tensor([0, 1, 2, 0, 3, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3])
    permutation = torch.randperm(values.shape[0], generator=generator)

    original, _ = segment_max(values, index, n_segments=4)
    permuted, _ = segment_max(values[permutation], index[permutation], n_segments=4)

    assert torch.allclose(original, permuted)


def test_segment_softmax_normalizes_each_group() -> None:
    logits = torch.tensor([1.0, 2.0, -1.0, 0.5, 0.5])
    index = torch.tensor([0, 0, 1, 2, 2])

    weights = segment_softmax(logits, index, n_segments=3)
    sums = segment_sum(weights, index, n_segments=3)

    assert torch.allclose(sums, torch.ones(3))
    assert weights[1] > weights[0]
    assert torch.allclose(weights[3:], torch.tensor([0.5, 0.5]))


def test_segment_sum_supports_gradients() -> None:
    values = torch.tensor([[1.0], [2.0], [3.0]], requires_grad=True)
    index = torch.tensor([0, 0, 1])

    output = segment_sum(values, index, n_segments=2)
    loss = output.square().sum()
    loss.backward()

    assert values.grad is not None
    assert torch.isfinite(values.grad).all()
