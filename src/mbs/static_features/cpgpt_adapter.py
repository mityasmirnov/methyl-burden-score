"""CpGPT sequence-adapter loader that avoids the broken torchtune/torchao import chain.

Reimplements the ``MLPBlock`` used as ``CpGPT.dna_encoder`` / ``encode_sequence()``
for the small (CpGPT2M) checkpoint so offline export works with MBS torch pins.
"""

from __future__ import annotations

import pickle
from pathlib import Path

import torch
import torch.nn.functional as F
from torch import nn


class _SwiGLU(nn.Module):
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x1, x2 = x.chunk(2, dim=-1)
        return x1 * F.silu(x2)


class SequenceAdapterMLP(nn.Module):
    """CpGPT ``MLPBlock`` used as the DNA sequence adapter (encode_sequence)."""

    def __init__(
        self,
        d_in: int,
        d_hidden: int,
        d_out: int,
        *,
        dropout: float,
        n_blocks: int,
        expansion_factor: int = 2,
        input_bias: bool = True,
        bias: bool = False,
        out_bias: bool = False,
        activation: str = "swiglu",
        pre_norm: bool = False,
        post_norm: bool = False,
        norm_type: str = "rmsnorm",
    ) -> None:
        super().__init__()
        if activation != "swiglu":
            raise ValueError(f"unsupported activation for Stage 0 export: {activation}")
        if norm_type != "rmsnorm":
            raise ValueError(f"unsupported norm_type for Stage 0 export: {norm_type}")

        self.input_norm: nn.Module = nn.RMSNorm(d_in) if pre_norm else nn.Identity()
        self.input_adapter = nn.Linear(d_in, d_hidden, bias=input_bias)
        self.blocks = nn.ModuleList()
        for _ in range(n_blocks):
            self.blocks.append(
                nn.Sequential(
                    nn.RMSNorm(d_hidden),
                    nn.Linear(
                        d_hidden,
                        expansion_factor * d_hidden * 2,
                        bias=bias,
                    ),
                    _SwiGLU(),
                    nn.Dropout(dropout),
                    nn.Linear(expansion_factor * d_hidden, d_hidden, bias=bias),
                )
            )
        self.output_norm: nn.Module = nn.RMSNorm(d_hidden) if post_norm else nn.Identity()
        self.output_adapter = nn.Linear(d_hidden, d_out, bias=out_bias)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = self.input_norm(x)
        x = self.input_adapter(x)
        for block in self.blocks:
            x = x + block(x)
        x = self.output_norm(x)
        return self.output_adapter(x)


class _Ignore:
    def __init__(self, *args: object, **kwargs: object) -> None:
        pass

    def __setstate__(self, state: object) -> None:
        pass


class _TolerantUnpickler(pickle.Unpickler):
    def find_class(self, module: str, name: str) -> object:
        try:
            return super().find_class(module, name)
        except Exception:
            return _Ignore


def load_checkpoint_state_dict(checkpoint_path: Path) -> dict[str, torch.Tensor]:
    """Load Lightning checkpoint ``state_dict`` without importing CpGPT modules."""
    path = checkpoint_path.resolve()
    if not path.is_file():
        raise FileNotFoundError(f"checkpoint not found: {path}")

    class _CompatPickle:
        Unpickler = _TolerantUnpickler

        @staticmethod
        def load(file: object, *args: object, **kwargs: object) -> object:
            return _TolerantUnpickler(file).load()  # type: ignore[arg-type]

    obj = torch.load(
        path,
        map_location="cpu",
        weights_only=False,
        pickle_module=_CompatPickle,
    )
    if not isinstance(obj, dict) or "state_dict" not in obj:
        raise ValueError(f"unexpected checkpoint payload in {path}")
    state_dict = obj["state_dict"]
    if not isinstance(state_dict, dict):
        raise TypeError(f"checkpoint state_dict is not a mapping: {path}")
    return state_dict


def load_small_sequence_adapter(
    checkpoint_path: Path,
    *,
    device: str = "cpu",
    dropout: float = 0.01,
    n_blocks: int = 3,
    d_in: int = 1024,
    d_hidden: int = 128,
    d_out: int = 128,
) -> SequenceAdapterMLP:
    """Instantiate and load the CpGPT2M (small) DNA sequence adapter."""
    state_dict = load_checkpoint_state_dict(checkpoint_path)
    mapped: dict[str, torch.Tensor] = {}
    for key, value in state_dict.items():
        key_s = str(key)
        if "dna_encoder." not in key_s:
            continue
        local_key = key_s.split("dna_encoder.", 1)[1]
        mapped[local_key] = value

    if not mapped:
        raise ValueError(f"no dna_encoder weights found in {checkpoint_path}")

    adapter = SequenceAdapterMLP(
        d_in,
        d_hidden,
        d_out,
        dropout=dropout,
        n_blocks=n_blocks,
        activation="swiglu",
        pre_norm=False,
        post_norm=False,
        bias=False,
        input_bias=True,
        out_bias=False,
        norm_type="rmsnorm",
    )
    missing, unexpected = adapter.load_state_dict(mapped, strict=True)
    if missing or unexpected:
        raise RuntimeError(f"dna_encoder load mismatch missing={missing} unexpected={unexpected}")
    adapter.to(torch.device(device))
    adapter.eval()
    return adapter
