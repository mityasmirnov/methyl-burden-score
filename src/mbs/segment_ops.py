"""Pure-PyTorch segment operations for ragged methylation sets."""

from __future__ import annotations

from typing import Literal

import torch
from torch import Tensor

PoolName = Literal["sum", "mean", "sqrt_sum", "max"]


def _validate(values: Tensor, segment_index: Tensor, n_segments: int) -> None:
    if values.ndim < 1:
        raise ValueError("values must have at least one dimension")
    if segment_index.ndim != 1:
        raise ValueError("segment_index must be one-dimensional")
    if values.shape[0] != segment_index.shape[0]:
        raise ValueError(
            "values and segment_index must have equal leading dimensions: "
            f"{values.shape[0]} != {segment_index.shape[0]}"
        )
    if n_segments < 0:
        raise ValueError("n_segments must be non-negative")
    if segment_index.dtype not in (torch.int32, torch.int64):
        raise TypeError("segment_index must have an integer dtype")
    if segment_index.numel() == 0:
        return
    minimum = int(segment_index.min().item())
    maximum = int(segment_index.max().item())
    if minimum < 0 or maximum >= n_segments:
        raise IndexError(
            f"segment indices must lie in [0, {n_segments}); observed [{minimum}, {maximum}]"
        )


def _expanded_index(values: Tensor, segment_index: Tensor) -> Tensor:
    shape = (segment_index.shape[0],) + (1,) * (values.ndim - 1)
    return segment_index.view(shape).expand_as(values)


def segment_count(segment_index: Tensor, n_segments: int) -> Tensor:
    """Count elements in each segment."""
    dummy = torch.empty(
        segment_index.shape[0],
        device=segment_index.device,
    )
    _validate(dummy, segment_index, n_segments)
    counts = torch.zeros(
        n_segments,
        dtype=torch.long,
        device=segment_index.device,
    )
    counts.index_add_(
        0,
        segment_index.to(torch.long),
        torch.ones_like(segment_index, dtype=torch.long),
    )
    return counts


def segment_sum(values: Tensor, segment_index: Tensor, n_segments: int) -> Tensor:
    """Sum values within segments."""
    _validate(values, segment_index, n_segments)
    output = torch.zeros(
        (n_segments, *values.shape[1:]),
        dtype=values.dtype,
        device=values.device,
    )
    output.index_add_(0, segment_index.to(torch.long), values)
    return output


def segment_mean(values: Tensor, segment_index: Tensor, n_segments: int) -> Tensor:
    """Mean values within segments, returning zero for empty segments."""
    total = segment_sum(values, segment_index, n_segments)
    counts = segment_count(segment_index, n_segments).to(values.dtype)
    denominator_shape = (n_segments,) + (1,) * (values.ndim - 1)
    denominator = counts.view(denominator_shape).clamp_min(1)
    return total / denominator


def segment_sqrt_sum(values: Tensor, segment_index: Tensor, n_segments: int) -> Tensor:
    """Sum divided by the square root of segment size."""
    total = segment_sum(values, segment_index, n_segments)
    counts = segment_count(segment_index, n_segments).to(values.dtype)
    denominator_shape = (n_segments,) + (1,) * (values.ndim - 1)
    denominator = counts.sqrt().view(denominator_shape).clamp_min(1)
    return total / denominator


def segment_max(
    values: Tensor,
    segment_index: Tensor,
    n_segments: int,
) -> tuple[Tensor, Tensor]:
    """Elementwise maximum and a presence mask for each segment."""
    _validate(values, segment_index, n_segments)
    if not values.is_floating_point():
        raise TypeError("segment_max requires floating-point values")

    output = torch.full(
        (n_segments, *values.shape[1:]),
        -torch.inf,
        dtype=values.dtype,
        device=values.device,
    )
    if values.shape[0] > 0:
        output.scatter_reduce_(
            0,
            _expanded_index(values, segment_index.to(torch.long)),
            values,
            reduce="amax",
            include_self=True,
        )

    present = segment_count(segment_index, n_segments) > 0
    mask_shape = (n_segments,) + (1,) * (values.ndim - 1)
    output = torch.where(
        present.view(mask_shape),
        output,
        torch.zeros_like(output),
    )
    return output, present


def segment_softmax(logits: Tensor, segment_index: Tensor, n_segments: int) -> Tensor:
    """Apply a numerically stable softmax independently within each segment."""
    if logits.ndim != 1:
        raise ValueError("segment_softmax currently expects one-dimensional logits")
    maxima, _ = segment_max(logits, segment_index, n_segments)
    shifted = logits - maxima[segment_index.to(torch.long)]
    numerator = shifted.exp()
    denominator = segment_sum(numerator, segment_index, n_segments)
    return numerator / denominator[segment_index.to(torch.long)].clamp_min(1e-12)


def segment_pool(
    values: Tensor,
    segment_index: Tensor,
    n_segments: int,
    pool: PoolName,
) -> tuple[Tensor, Tensor]:
    """Pool values and return a segment-presence mask."""
    present = segment_count(segment_index, n_segments) > 0
    if pool == "sum":
        return segment_sum(values, segment_index, n_segments), present
    if pool == "mean":
        return segment_mean(values, segment_index, n_segments), present
    if pool == "sqrt_sum":
        return segment_sqrt_sum(values, segment_index, n_segments), present
    if pool == "max":
        return segment_max(values, segment_index, n_segments)
    raise ValueError(f"unknown pooling operator: {pool}")
