"""Stage 0 permutation-invariant methylation burden models."""

from __future__ import annotations

from collections.abc import Sequence
from typing import TypedDict

import torch
import torch.nn.functional as F
from torch import Tensor, nn

from mbs.segment_ops import PoolName, segment_pool


class ModelOutput(TypedDict):
    mbs: Tensor
    centered_mbs: Tensor
    present: Tensor
    logits: Tensor


class HierarchicalModelOutput(ModelOutput):
    region_hidden: Tensor
    region_present: Tensor
    gene_hidden: Tensor


class SharedMLP(nn.Module):
    """Small MLP used as a shared element or set encoder."""

    def __init__(
        self,
        input_dim: int,
        hidden_dims: Sequence[int],
        output_dim: int,
        *,
        dropout: float = 0.0,
        layer_norm: bool = False,
        activation: str = "leaky_relu",
    ) -> None:
        super().__init__()
        if input_dim <= 0 or output_dim <= 0:
            raise ValueError("input_dim and output_dim must be positive")

        dimensions = [input_dim, *hidden_dims, output_dim]
        layers: list[nn.Module] = []
        for index, (source, target) in enumerate(zip(dimensions[:-1], dimensions[1:], strict=True)):
            layers.append(nn.Linear(source, target))
            is_last = index == len(dimensions) - 2
            if is_last:
                continue
            if layer_norm:
                layers.append(nn.LayerNorm(target))
            if activation == "leaky_relu":
                layers.append(nn.LeakyReLU(negative_slope=0.01))
            elif activation == "gelu":
                layers.append(nn.GELU())
            else:
                raise ValueError(f"unsupported activation: {activation}")
            if dropout > 0:
                layers.append(nn.Dropout(dropout))
        self.network = nn.Sequential(*layers)

    def forward(self, inputs: Tensor) -> Tensor:
        return self.network(inputs)


class FlatDeepSet(nn.Module):
    """Direct CpG-to-gene Deep Set reference model."""

    def __init__(
        self,
        input_dim: int,
        *,
        phi_hidden_dim: int = 20,
        phi_layers: int = 2,
        rho_hidden_dim: int = 10,
        rho_layers: int = 3,
        pool: PoolName = "max",
        neutral_score: float = 0.5,
        dropout: float = 0.0,
    ) -> None:
        super().__init__()
        if phi_layers < 1 or rho_layers < 1:
            raise ValueError("phi_layers and rho_layers must be at least one")
        self.pool = pool
        self.neutral_score = neutral_score
        self.phi = SharedMLP(
            input_dim,
            [phi_hidden_dim] * (phi_layers - 1),
            phi_hidden_dim,
            dropout=dropout,
            activation="leaky_relu",
        )
        self.rho = SharedMLP(
            phi_hidden_dim,
            [rho_hidden_dim] * (rho_layers - 1),
            1,
            dropout=dropout,
            activation="leaky_relu",
        )

    def forward(
        self,
        cpg_features: Tensor,
        cpg_to_gene: Tensor,
        n_gene_instances: int,
    ) -> ModelOutput:
        cpg_hidden = self.phi(cpg_features)
        gene_hidden, present = segment_pool(
            cpg_hidden,
            cpg_to_gene,
            n_gene_instances,
            self.pool,
        )
        logits = self.rho(gene_hidden).squeeze(-1)
        raw_score = torch.sigmoid(logits)
        neutral = torch.full_like(raw_score, self.neutral_score)
        mbs = torch.where(present, raw_score, neutral)
        centered = torch.where(
            present,
            mbs - self.neutral_score,
            torch.zeros_like(mbs),
        )
        return {
            "mbs": mbs,
            "centered_mbs": centered,
            "present": present,
            "logits": logits,
        }


class HierarchicalDeepSet(nn.Module):
    """Shared CpG-to-region-to-gene methylation burden model."""

    def __init__(
        self,
        input_dim: int,
        n_region_types: int,
        *,
        cpg_hidden_dim: int = 64,
        region_hidden_dim: int = 32,
        region_type_dim: int = 8,
        cpg_pool: PoolName = "max",
        region_pool: PoolName = "max",
        neutral_score: float = 0.5,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        if n_region_types <= 0:
            raise ValueError("n_region_types must be positive")
        self.cpg_pool = cpg_pool
        self.region_pool = region_pool
        self.neutral_score = neutral_score

        self.cpg_encoder = SharedMLP(
            input_dim,
            [cpg_hidden_dim],
            cpg_hidden_dim,
            dropout=dropout,
            layer_norm=True,
            activation="gelu",
        )
        self.region_type_embedding = nn.Embedding(n_region_types, region_type_dim)
        self.region_encoder = SharedMLP(
            cpg_hidden_dim + region_type_dim,
            [region_hidden_dim],
            region_hidden_dim,
            dropout=dropout,
            layer_norm=True,
            activation="gelu",
        )
        self.rho = SharedMLP(
            region_hidden_dim,
            [16, 8],
            1,
            dropout=dropout,
            activation="leaky_relu",
        )

    def forward(
        self,
        *,
        cpg_features: Tensor,
        cpg_to_region: Tensor,
        region_type: Tensor,
        region_to_gene: Tensor,
        n_regions: int,
        n_gene_instances: int,
    ) -> HierarchicalModelOutput:
        if region_type.shape != (n_regions,):
            raise ValueError(
                f"region_type must have shape ({n_regions},), found {tuple(region_type.shape)}"
            )
        if region_to_gene.shape != (n_regions,):
            raise ValueError(
                "region_to_gene must have one entry per region: "
                f"expected {n_regions}, found {region_to_gene.shape[0]}"
            )

        cpg_hidden = self.cpg_encoder(cpg_features)
        region_pooled, region_present = segment_pool(
            cpg_hidden,
            cpg_to_region,
            n_regions,
            self.cpg_pool,
        )
        type_hidden = self.region_type_embedding(region_type.to(torch.long))
        region_hidden = self.region_encoder(torch.cat([region_pooled, type_hidden], dim=-1))
        region_hidden = region_hidden * region_present.unsqueeze(-1)

        active_region_hidden = region_hidden[region_present]
        active_region_to_gene = region_to_gene[region_present]
        gene_hidden, gene_present = segment_pool(
            active_region_hidden,
            active_region_to_gene,
            n_gene_instances,
            self.region_pool,
        )

        logits = self.rho(gene_hidden).squeeze(-1)
        raw_score = torch.sigmoid(logits)
        neutral = torch.full_like(raw_score, self.neutral_score)
        mbs = torch.where(gene_present, raw_score, neutral)
        centered = torch.where(
            gene_present,
            mbs - self.neutral_score,
            torch.zeros_like(mbs),
        )

        return {
            "mbs": mbs,
            "centered_mbs": centered,
            "present": gene_present,
            "logits": logits,
            "region_hidden": region_hidden,
            "region_present": region_present,
            "gene_hidden": gene_hidden,
        }


class SeedMaskedLinearHead(nn.Module):
    """Trait-specific linear head over shared gene scores and optional covariates."""

    def __init__(
        self,
        n_genes: int,
        n_outputs: int,
        seed_mask: Tensor,
        *,
        n_covariates: int = 0,
        neutral_score: float = 0.5,
    ) -> None:
        super().__init__()
        expected_shape = (n_outputs, n_genes)
        if tuple(seed_mask.shape) != expected_shape:
            raise ValueError(
                f"seed_mask must have shape {expected_shape}, found {tuple(seed_mask.shape)}"
            )
        self.n_genes = n_genes
        self.n_covariates = n_covariates
        self.neutral_score = neutral_score
        self.register_buffer("seed_mask", seed_mask.to(torch.float32))
        self.gene_weight = nn.Parameter(torch.zeros(n_outputs, n_genes))
        self.bias = nn.Parameter(torch.zeros(n_outputs))
        self.covariate_head = (
            nn.Linear(n_covariates, n_outputs, bias=False)
            if n_covariates > 0
            else None
        )

    def forward(
        self,
        mbs: Tensor,
        present: Tensor,
        covariates: Tensor | None = None,
    ) -> Tensor:
        if mbs.ndim != 2 or mbs.shape[1] != self.n_genes:
            raise ValueError(f"mbs must have shape [batch, {self.n_genes}]")
        if present.shape != mbs.shape:
            raise ValueError("present mask must match mbs shape")

        centered = torch.where(
            present,
            mbs - self.neutral_score,
            torch.zeros_like(mbs),
        )
        output = F.linear(
            centered,
            self.gene_weight * self.seed_mask,
            self.bias,
        )

        if self.covariate_head is not None:
            if covariates is None:
                raise ValueError("covariates are required for this head")
            if covariates.ndim != 2 or covariates.shape[1] != self.n_covariates:
                raise ValueError(
                    f"covariates must have shape [batch, {self.n_covariates}]"
                )
            output = output + self.covariate_head(covariates)
        elif covariates is not None:
            raise ValueError("covariates were supplied to a head configured without them")

        return output
