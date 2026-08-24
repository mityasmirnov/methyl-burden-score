"""Shared encoder hyperparameters for parameter-matched flat vs hierarchical."""

from __future__ import annotations

from typing import Any


def resolve_encoder(
    model_cfg: dict[str, Any],
    *,
    default_activation: str = "leaky_relu",
    default_layer_norm: bool = False,
    default_dropout: float = 0.0,
    default_cpg_hidden: int = 20,
) -> dict[str, Any]:
    """Read ``model.encoder`` with fallbacks to legacy flat/hier keys."""
    raw_enc = model_cfg.get("encoder")
    enc: dict[str, Any] = raw_enc if isinstance(raw_enc, dict) else {}
    dropout = float(enc.get("dropout", model_cfg.get("dropout", default_dropout)))
    activation = str(enc.get("activation", model_cfg.get("activation", default_activation)))
    layer_norm = bool(enc.get("layer_norm", model_cfg.get("layer_norm", default_layer_norm)))
    cpg_hidden = int(
        enc.get(
            "cpg_hidden_dim",
            model_cfg.get(
                "cpg_hidden_dimension",
                model_cfg.get("phi_hidden_dimension", default_cpg_hidden),
            ),
        )
    )
    return {
        "activation": activation,
        "dropout": dropout,
        "layer_norm": layer_norm,
        "cpg_hidden_dim": cpg_hidden,
    }
