"""Ragged batch contract for Stage 0 methylation models."""

from __future__ import annotations

from dataclasses import dataclass, fields

import torch
from torch import Tensor


@dataclass(slots=True)
class MethylationBatch:
    """Observed sample–CpG pairs and their hierarchical segment mappings."""

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

    def validate(self) -> None:
        """Validate dimensional and index contracts without changing the batch."""
        n_cpg = self.cpg_features.shape[0]
        for name, tensor in (
            ("locus_row", self.locus_row),
            ("cpg_sample_index", self.cpg_sample_index),
            ("cpg_to_region", self.cpg_to_region),
        ):
            if tensor.ndim != 1 or tensor.shape[0] != n_cpg:
                raise ValueError(f"{name} must have shape [{n_cpg}]")

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
