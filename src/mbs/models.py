"""Stage 0 permutation-invariant methylation burden models."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Literal, TypedDict

import torch
import torch.nn.functional as F
from torch import Tensor, nn

from mbs.segment_ops import PoolName, segment_pool


def center_mask_scores(
    mbs: Tensor,
    present: Tensor,
    *,
    neutral_score: float = 0.5,
) -> Tensor:
    """Absent scores contribute 0 to a linear head; present scores are centered."""
    return torch.where(
        present.to(dtype=torch.bool),
        mbs - neutral_score,
        torch.zeros_like(mbs),
    )


class ModelOutput(TypedDict):
    mbs: Tensor
    centered_mbs: Tensor
    present: Tensor
    logits: Tensor


class HierarchicalModelOutput(ModelOutput):
    region_hidden: Tensor
    region_present: Tensor
    gene_hidden: Tensor
    residual_mbs: Tensor
    residual_present: Tensor
    residual_logits: Tensor


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
        activation: str = "leaky_relu",
        layer_norm: bool = False,
    ) -> None:
        super().__init__()
        if phi_layers < 1 or rho_layers < 1:
            raise ValueError("phi_layers and rho_layers must be at least one")
        self.pool: PoolName = pool
        self.neutral_score = neutral_score
        self.phi = SharedMLP(
            input_dim,
            [phi_hidden_dim] * (phi_layers - 1),
            phi_hidden_dim,
            dropout=dropout,
            layer_norm=layer_norm,
            activation=activation,
        )
        self.rho = SharedMLP(
            phi_hidden_dim,
            [rho_hidden_dim] * (rho_layers - 1),
            1,
            dropout=dropout,
            layer_norm=layer_norm,
            activation=activation,
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
        centered = center_mask_scores(mbs, present, neutral_score=self.neutral_score)
        return {
            "mbs": mbs,
            "centered_mbs": centered,
            "present": present,
            "logits": logits,
        }


class FlatDeepSetRegion(FlatDeepSet):
    """Annotation-augmented flat pooling: per-CpG features → gene MBS.

    Programme names: Stage A **N-light-gene-max** / **N-light-gene-mean**.
    Input layout is assembled in ``training.flat_region_features``
    (M-value, gene-role, CGI context, regulatory multi-hot, presence flags,
    observed) — no RBS intermediate hop.
    """


class HierarchicalDeepSet(nn.Module):
    """CpG→region→gene hierarchy for annotated loci + residual path for unmapped.

    Unmapped / ambiguous loci are **not** pooled into genes. They use a separate
    residual DeepSet (shared CpG encoder → per-sample max pool → residual score).
    """

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
        residual_pool: PoolName = "max",
        neutral_score: float = 0.5,
        dropout: float = 0.1,
        activation: str = "gelu",
        layer_norm: bool = True,
    ) -> None:
        super().__init__()
        if n_region_types <= 0:
            raise ValueError("n_region_types must be positive")
        self.cpg_pool: PoolName = cpg_pool
        self.region_pool: PoolName = region_pool
        self.residual_pool: PoolName = residual_pool
        self.neutral_score = neutral_score
        self.cpg_hidden_dim = cpg_hidden_dim
        self.region_hidden_dim = region_hidden_dim

        self.cpg_encoder = SharedMLP(
            input_dim,
            [cpg_hidden_dim],
            cpg_hidden_dim,
            dropout=dropout,
            layer_norm=layer_norm,
            activation=activation,
        )
        self.region_type_embedding = nn.Embedding(n_region_types, region_type_dim)
        self.region_encoder = SharedMLP(
            cpg_hidden_dim + region_type_dim,
            [region_hidden_dim],
            region_hidden_dim,
            dropout=dropout,
            layer_norm=layer_norm,
            activation=activation,
        )
        self.rho = SharedMLP(
            region_hidden_dim,
            [16, 8],
            1,
            dropout=dropout,
            activation=activation,
        )
        self.residual_rho = SharedMLP(
            cpg_hidden_dim,
            [16, 8],
            1,
            dropout=dropout,
            activation=activation,
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
        residual_features: Tensor | None = None,
        residual_sample_index: Tensor | None = None,
        n_samples: int | None = None,
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

        device = cpg_features.device if cpg_features.numel() else region_type.device
        if cpg_features.shape[0] == 0:
            region_hidden = torch.zeros(
                n_regions, self.region_hidden_dim, device=device, dtype=torch.float32
            )
            region_present = torch.zeros(n_regions, dtype=torch.bool, device=device)
            gene_hidden = torch.zeros(
                n_gene_instances, self.region_hidden_dim, device=device, dtype=torch.float32
            )
            gene_present = torch.zeros(n_gene_instances, dtype=torch.bool, device=device)
            logits = torch.zeros(n_gene_instances, device=device, dtype=torch.float32)
            mbs = torch.full_like(logits, self.neutral_score)
            centered = torch.zeros_like(mbs)
        else:
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
            centered = center_mask_scores(mbs, gene_present, neutral_score=self.neutral_score)

        residual_mbs, residual_present, residual_logits = self._forward_residual(
            residual_features=residual_features,
            residual_sample_index=residual_sample_index,
            n_samples=n_samples,
            device=mbs.device,
            dtype=mbs.dtype,
        )

        return {
            "mbs": mbs,
            "centered_mbs": centered,
            "present": gene_present,
            "logits": logits,
            "region_hidden": region_hidden,
            "region_present": region_present,
            "gene_hidden": gene_hidden,
            "residual_mbs": residual_mbs,
            "residual_present": residual_present,
            "residual_logits": residual_logits,
        }

    def _forward_residual(
        self,
        *,
        residual_features: Tensor | None,
        residual_sample_index: Tensor | None,
        n_samples: int | None,
        device: torch.device,
        dtype: torch.dtype,
    ) -> tuple[Tensor, Tensor, Tensor]:
        if n_samples is None:
            n_samples = 0
        if residual_features is None or residual_sample_index is None or n_samples <= 0:
            zeros = torch.zeros(n_samples, device=device, dtype=dtype)
            present = torch.zeros(n_samples, dtype=torch.bool, device=device)
            return zeros, present, zeros
        if residual_features.shape[0] != residual_sample_index.shape[0]:
            raise ValueError("residual_features and residual_sample_index length mismatch")
        if residual_features.shape[0] == 0:
            zeros = torch.zeros(n_samples, device=device, dtype=dtype)
            present = torch.zeros(n_samples, dtype=torch.bool, device=device)
            return zeros, present, zeros

        residual_hidden = self.cpg_encoder(residual_features)
        pooled, present = segment_pool(
            residual_hidden,
            residual_sample_index,
            n_samples,
            self.residual_pool,
        )
        logits = self.residual_rho(pooled).squeeze(-1)
        raw = torch.sigmoid(logits)
        neutral = torch.full_like(raw, self.neutral_score)
        mbs = torch.where(present, raw, neutral)
        return mbs, present, logits


class CascadeModelOutput(TypedDict):
    rbs: Tensor
    rbs_present: Tensor
    region_hidden: Tensor
    mbs: Tensor
    centered_mbs: Tensor
    present: Tensor
    orphan_rbs: Tensor
    orphan_present: Tensor
    logits: Tensor


GeneAggregation = Literal["scalar_rbs", "region_hidden"]


class CascadeDeepSet(nn.Module):
    """7F topology: CpG → typed region → gene-pooled MBS; orphan RBS kept.

    No residual/tile path — leftover CpGs are scored outside this module
    (direct elastic-net). ``region_to_gene`` uses -1 for orphan regions.

    ``gene_aggregation``:
    - ``scalar_rbs``: pool sigmoid scalar RBS by gene (legacy product path).
    - ``region_hidden``: pool latent region embeddings by gene, then ``gene_rho``
      → MBS (DeepRVAT-like; scalar RBS still exported for diagnostics).
    """

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
        gene_aggregation: GeneAggregation = "scalar_rbs",
        neutral_score: float = 0.5,
        dropout: float = 0.1,
        activation: str = "gelu",
        layer_norm: bool = True,
    ) -> None:
        super().__init__()
        if n_region_types <= 0:
            raise ValueError("n_region_types must be positive")
        if gene_aggregation not in ("scalar_rbs", "region_hidden"):
            raise ValueError(f"unsupported gene_aggregation: {gene_aggregation!r}")
        self.cpg_pool: PoolName = cpg_pool
        self.region_pool: PoolName = region_pool
        self.gene_aggregation: GeneAggregation = gene_aggregation
        self.neutral_score = neutral_score
        self.cpg_hidden_dim = cpg_hidden_dim
        self.region_hidden_dim = region_hidden_dim

        self.cpg_encoder = SharedMLP(
            input_dim,
            [cpg_hidden_dim],
            cpg_hidden_dim,
            dropout=dropout,
            layer_norm=layer_norm,
            activation=activation,
        )
        self.region_type_embedding = nn.Embedding(n_region_types, region_type_dim)
        self.region_encoder = SharedMLP(
            cpg_hidden_dim + region_type_dim,
            [region_hidden_dim],
            region_hidden_dim,
            dropout=dropout,
            layer_norm=layer_norm,
            activation=activation,
        )
        self.region_rho = SharedMLP(
            region_hidden_dim,
            [region_hidden_dim],
            1,
            dropout=dropout,
            activation=activation,
        )
        # Only for region_hidden; keep off scalar_rbs so old P2/P4 checkpoints load.
        self.gene_rho: SharedMLP | None
        if gene_aggregation == "region_hidden":
            self.gene_rho = SharedMLP(
                region_hidden_dim,
                [16, 8],
                1,
                dropout=dropout,
                activation=activation,
            )
        else:
            self.gene_rho = None

    def forward_from_cpg_hidden(
        self,
        cpg_hidden: Tensor,
        *,
        cpg_to_region: Tensor,
        region_type: Tensor,
        region_to_gene: Tensor,
        n_regions: int,
        n_gene_instances: int,
        orphan_region_indices: Tensor | None = None,
    ) -> CascadeModelOutput:
        """Region→gene path given per-edge CpG encoder outputs (one sample)."""
        if region_type.shape != (n_regions,):
            raise ValueError(
                f"region_type must have shape ({n_regions},), found {tuple(region_type.shape)}"
            )
        if region_to_gene.shape != (n_regions,):
            raise ValueError(
                "region_to_gene must have one entry per region: "
                f"expected {n_regions}, found {region_to_gene.shape[0]}"
            )

        device = cpg_hidden.device if cpg_hidden.numel() else region_type.device
        dtype = cpg_hidden.dtype if cpg_hidden.numel() else torch.float32

        if cpg_hidden.shape[0] == 0 or n_regions == 0:
            rbs = torch.full((n_regions,), self.neutral_score, device=device, dtype=dtype)
            rbs_present = torch.zeros(n_regions, dtype=torch.bool, device=device)
            region_hidden = torch.zeros(
                n_regions, self.region_hidden_dim, device=device, dtype=dtype
            )
            mbs = torch.full((n_gene_instances,), self.neutral_score, device=device, dtype=dtype)
            present = torch.zeros(n_gene_instances, dtype=torch.bool, device=device)
            logits = torch.zeros(n_gene_instances, device=device, dtype=dtype)
        else:
            region_pooled, rbs_present = segment_pool(
                cpg_hidden,
                cpg_to_region,
                n_regions,
                self.cpg_pool,
            )
            type_hidden = self.region_type_embedding(region_type.to(torch.long))
            region_hidden = self.region_encoder(torch.cat([region_pooled, type_hidden], dim=-1))
            region_hidden = region_hidden * rbs_present.unsqueeze(-1)
            rbs_logits = self.region_rho(region_hidden).squeeze(-1)
            raw_rbs = torch.sigmoid(rbs_logits)
            rbs = torch.where(
                rbs_present,
                raw_rbs,
                torch.full_like(raw_rbs, self.neutral_score),
            )

            allocated = region_to_gene >= 0
            if allocated.any() and n_gene_instances > 0:
                active = allocated & rbs_present
                if active.any():
                    gene_idx = region_to_gene[active].to(torch.long)
                    if self.gene_aggregation == "region_hidden":
                        if self.gene_rho is None:
                            raise RuntimeError("gene_rho required for region_hidden aggregation")
                        gene_hidden, gene_present = segment_pool(
                            region_hidden[active],
                            gene_idx,
                            n_gene_instances,
                            self.region_pool,
                        )
                        logits = self.gene_rho(gene_hidden).squeeze(-1)
                        raw_mbs = torch.sigmoid(logits)
                        mbs = torch.where(
                            gene_present,
                            raw_mbs,
                            torch.full_like(raw_mbs, self.neutral_score),
                        )
                        present = gene_present
                    else:
                        score_vals = rbs[active].unsqueeze(-1)
                        gene_scores, gene_present = segment_pool(
                            score_vals,
                            gene_idx,
                            n_gene_instances,
                            self.region_pool,
                        )
                        mbs = gene_scores.squeeze(-1)
                        mbs = torch.where(
                            gene_present,
                            mbs,
                            torch.full_like(mbs, self.neutral_score),
                        )
                        present = gene_present
                        logits = mbs
                else:
                    mbs = torch.full(
                        (n_gene_instances,), self.neutral_score, device=device, dtype=dtype
                    )
                    present = torch.zeros(n_gene_instances, dtype=torch.bool, device=device)
                    logits = torch.zeros(n_gene_instances, device=device, dtype=dtype)
            else:
                mbs = torch.full(
                    (n_gene_instances,), self.neutral_score, device=device, dtype=dtype
                )
                present = torch.zeros(n_gene_instances, dtype=torch.bool, device=device)
                logits = torch.zeros(n_gene_instances, device=device, dtype=dtype)

        if orphan_region_indices is None:
            orphan_idx = torch.nonzero(region_to_gene < 0, as_tuple=False).flatten()
        else:
            orphan_idx = orphan_region_indices.to(torch.long)
        if orphan_idx.numel() == 0:
            orphan_rbs = torch.zeros(0, device=device, dtype=dtype)
            orphan_present = torch.zeros(0, dtype=torch.bool, device=device)
        else:
            orphan_rbs = rbs[orphan_idx]
            orphan_present = rbs_present[orphan_idx]

        centered = center_mask_scores(mbs, present, neutral_score=self.neutral_score)
        return {
            "rbs": rbs,
            "rbs_present": rbs_present,
            "region_hidden": region_hidden,
            "mbs": mbs,
            "centered_mbs": centered,
            "present": present,
            "orphan_rbs": orphan_rbs,
            "orphan_present": orphan_present,
            "logits": logits,
        }

    def forward(
        self,
        *,
        cpg_features: Tensor,
        cpg_to_region: Tensor,
        region_type: Tensor,
        region_to_gene: Tensor,
        n_regions: int,
        n_gene_instances: int,
        orphan_region_indices: Tensor | None = None,
    ) -> CascadeModelOutput:
        if region_type.shape != (n_regions,):
            raise ValueError(
                f"region_type must have shape ({n_regions},), found {tuple(region_type.shape)}"
            )
        if region_to_gene.shape != (n_regions,):
            raise ValueError(
                "region_to_gene must have one entry per region: "
                f"expected {n_regions}, found {region_to_gene.shape[0]}"
            )

        device = cpg_features.device if cpg_features.numel() else region_type.device
        dtype = cpg_features.dtype if cpg_features.numel() else torch.float32

        if cpg_features.shape[0] == 0 or n_regions == 0:
            return self.forward_from_cpg_hidden(
                torch.zeros(0, self.cpg_hidden_dim, device=device, dtype=dtype),
                cpg_to_region=cpg_to_region,
                region_type=region_type,
                region_to_gene=region_to_gene,
                n_regions=n_regions,
                n_gene_instances=n_gene_instances,
                orphan_region_indices=orphan_region_indices,
            )

        cpg_hidden = self.cpg_encoder(cpg_features)
        return self.forward_from_cpg_hidden(
            cpg_hidden,
            cpg_to_region=cpg_to_region,
            region_type=region_type,
            region_to_gene=region_to_gene,
            n_regions=n_regions,
            n_gene_instances=n_gene_instances,
            orphan_region_indices=orphan_region_indices,
        )


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
            nn.Linear(n_covariates, n_outputs, bias=False) if n_covariates > 0 else None
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

        centered = center_mask_scores(mbs, present, neutral_score=self.neutral_score)
        seed_mask = self.get_buffer("seed_mask")
        output = F.linear(
            centered,
            self.gene_weight * seed_mask,
            self.bias,
        )

        if self.covariate_head is not None:
            if covariates is None:
                raise ValueError("covariates are required for this head")
            if covariates.ndim != 2 or covariates.shape[1] != self.n_covariates:
                raise ValueError(f"covariates must have shape [batch, {self.n_covariates}]")
            output = output + self.covariate_head(covariates)
        elif covariates is not None:
            raise ValueError("covariates were supplied to a head configured without them")

        return output
