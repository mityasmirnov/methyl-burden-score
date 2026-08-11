"""Ragged batch contract for Stage 0 methylation models."""

from __future__ import annotations

from dataclasses import dataclass, fields
from typing import Final

import torch
from torch import Tensor

# Per observed CpG / matrix column annotation status (Milestone 6).
ANNOTATION_STATUS_MAPPED: Final[int] = 0
ANNOTATION_STATUS_UNMAPPED: Final[int] = 1
ANNOTATION_STATUS_AMBIGUOUS: Final[int] = 2
ANNOTATION_STATUS_MULTI_MAPPED: Final[int] = 3

ANNOTATION_STATUS_NAMES: Final[tuple[str, ...]] = (
    "mapped",
    "unmapped",
    "ambiguous",
    "multi_mapped",
)

ANNOTATION_STATUS_TO_ID: Final[dict[str, int]] = {
    name: i for i, name in enumerate(ANNOTATION_STATUS_NAMES)
}


def annotation_status_masks(status: Tensor) -> dict[str, Tensor]:
    """Boolean masks for each annotation status (same length as ``status``)."""
    if status.ndim != 1:
        raise ValueError("annotation status must be one-dimensional")
    status_long = status.to(torch.long)
    return {
        "mapped": status_long == ANNOTATION_STATUS_MAPPED,
        "unmapped": status_long == ANNOTATION_STATUS_UNMAPPED,
        "ambiguous": status_long == ANNOTATION_STATUS_AMBIGUOUS,
        "multi_mapped": status_long == ANNOTATION_STATUS_MULTI_MAPPED,
    }


@dataclass(slots=True)
class MethylationBatch:
    """Observed sample–CpG pairs and their hierarchical segment mappings.

    Annotation-status tensors are aligned with observed CpG rows. Residual
    (unmapped / ambiguous) CpGs use ``residual_*`` fields and must not enter
    ``cpg_to_region`` gene pooling.
    """

    sample_ids: list[str]
    cpg_features: Tensor
    locus_row: Tensor
    cpg_sample_index: Tensor
    cpg_to_region: Tensor
    region_type: Tensor
    region_to_gene: Tensor
    gene_sample_index: Tensor
    gene_panel_index: Tensor
    targets: dict[str, Tensor]
    target_masks: dict[str, Tensor]
    covariates: dict[str, Tensor]
    annotation_status: Tensor
    residual_features: Tensor | None = None
    residual_sample_index: Tensor | None = None

    def annotation_masks(self) -> dict[str, Tensor]:
        """Return mapped / unmapped / ambiguous / multi_mapped masks."""
        return annotation_status_masks(self.annotation_status)

    def validate(self) -> None:
        """Validate dimensional and index contracts without changing the batch."""
        n_cpg = self.cpg_features.shape[0]
        for name, tensor in (
            ("locus_row", self.locus_row),
            ("cpg_sample_index", self.cpg_sample_index),
            ("cpg_to_region", self.cpg_to_region),
            ("annotation_status", self.annotation_status),
        ):
            if tensor.ndim != 1 or tensor.shape[0] != n_cpg:
                raise ValueError(f"{name} must have shape [{n_cpg}]")

        if self.annotation_status.numel() > 0:
            status_min = int(self.annotation_status.min().item())
            status_max = int(self.annotation_status.max().item())
            if status_min < 0 or status_max >= len(ANNOTATION_STATUS_NAMES):
                raise ValueError(
                    f"annotation_status values must be in [0, {len(ANNOTATION_STATUS_NAMES) - 1}]"
                )

        n_regions = self.region_type.shape[0]
        if self.region_type.ndim != 1:
            raise ValueError("region_type must be one-dimensional")
        if self.region_to_gene.ndim != 1 or self.region_to_gene.shape[0] != n_regions:
            raise ValueError("region_to_gene must have one entry per region")

        n_gene_instances = self.gene_sample_index.shape[0]
        if self.gene_sample_index.ndim != 1:
            raise ValueError("gene_sample_index must be one-dimensional")
        if self.gene_panel_index.shape != (n_gene_instances,):
            raise ValueError("gene_panel_index must have one entry per gene instance")

        if self.cpg_sample_index.numel() > 0:
            maximum_sample = int(self.cpg_sample_index.max().item())
            if maximum_sample >= len(self.sample_ids):
                raise IndexError("cpg_sample_index references a missing sample")

        if self.cpg_to_region.numel() > 0:
            maximum_region = int(self.cpg_to_region.max().item())
            if maximum_region >= n_regions:
                raise IndexError("cpg_to_region references a missing region")

        if self.region_to_gene.numel() > 0:
            maximum_gene = int(self.region_to_gene.max().item())
            if maximum_gene >= n_gene_instances:
                raise IndexError("region_to_gene references a missing gene instance")

        residual_features = self.residual_features
        residual_sample_index = self.residual_sample_index
        if (residual_features is None) != (residual_sample_index is None):
            raise ValueError(
                "residual_features and residual_sample_index must both be set or both None"
            )
        if residual_features is not None and residual_sample_index is not None:
            if residual_features.ndim != 2:
                raise ValueError("residual_features must have shape [N_residual, D]")
            if (
                residual_sample_index.ndim != 1
                or residual_sample_index.shape[0] != residual_features.shape[0]
            ):
                raise ValueError("residual_sample_index must have one entry per residual row")
            if residual_sample_index.numel() > 0:
                maximum_sample = int(residual_sample_index.max().item())
                if maximum_sample >= len(self.sample_ids):
                    raise IndexError("residual_sample_index references a missing sample")

        for task, target in self.targets.items():
            if task not in self.target_masks:
                raise KeyError(f"target mask is missing for task {task!r}")
            if self.target_masks[task].shape != target.shape:
                raise ValueError(f"target mask shape differs for task {task!r}")

    def to(self, device: torch.device | str, *, non_blocking: bool = False) -> MethylationBatch:
        """Move tensor fields to a device and preserve sample identifiers."""
        values: dict[str, object] = {}
        for field in fields(self):
            value = getattr(self, field.name)
            if isinstance(value, Tensor):
                values[field.name] = value.to(device, non_blocking=non_blocking)
            elif isinstance(value, dict):
                values[field.name] = {
                    key: tensor.to(device, non_blocking=non_blocking)
                    for key, tensor in value.items()
                }
            else:
                values[field.name] = value
        moved = MethylationBatch(**values)  # type: ignore[arg-type]
        moved.validate()
        return moved
