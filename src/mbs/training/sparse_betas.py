"""Minibatch row×column beta gathers with IO instrumentation (§3.5)."""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Any

import numpy as np


@dataclass
class GatherIoStats:
    """Accumulated gather metrics for one epoch / run."""

    n_calls: int = 0
    n_rows: int = 0
    n_logical_cols: int = 0
    n_unique_logical_cols: int = 0
    bytes_out: int = 0
    io_seconds: float = 0.0
    pack_seconds: float = 0.0

    def as_dict(self) -> dict[str, Any]:
        return {
            "n_calls": self.n_calls,
            "n_rows": self.n_rows,
            "n_logical_cols": self.n_logical_cols,
            "n_unique_logical_cols": self.n_unique_logical_cols,
            "bytes_out": self.bytes_out,
            "io_seconds": self.io_seconds,
            "pack_seconds": self.pack_seconds,
            "mean_io_ms": (1000.0 * self.io_seconds / self.n_calls) if self.n_calls else 0.0,
        }


@dataclass
class SparseBetasView:
    """Thin wrapper: fetch only requested columns; never preload the universe."""

    handle: Any
    n_cols: int
    stats: GatherIoStats = field(default_factory=GatherIoStats)

    def gather_rows_cols(
        self,
        row_ids: np.ndarray,
        col_ids: np.ndarray,
    ) -> np.ndarray:
        """Return ``(n_rows, n_cols_req)`` float32 block for ``col_ids`` order."""
        rows = np.asarray(row_ids, dtype=np.int64).reshape(-1)
        cols = np.asarray(col_ids, dtype=np.int64).reshape(-1)
        t0 = time.perf_counter()
        block = np.asarray(self.handle[np.ix_(rows, cols)], dtype=np.float32)
        io_s = time.perf_counter() - t0
        if block.shape != (rows.size, cols.size):
            raise RuntimeError(
                f"gather shape {block.shape} != expected {(rows.size, cols.size)}"
            )
        self.stats.n_calls += 1
        self.stats.n_rows += int(rows.size)
        self.stats.n_logical_cols += int(cols.size)
        self.stats.n_unique_logical_cols += int(np.unique(cols).size)
        self.stats.bytes_out += int(block.nbytes)
        self.stats.io_seconds += float(io_s)
        return block

    def beta_edge_values(self, row: int, edge_col_index: np.ndarray) -> np.ndarray:
        """Betas aligned to ``edge_col_index`` (unique fetch + restore)."""
        cols = np.asarray(edge_col_index, dtype=np.int64)
        t_pack0 = time.perf_counter()
        uniq, restore = np.unique(cols, return_inverse=True)
        self.stats.pack_seconds += time.perf_counter() - t_pack0
        fetched = self.gather_rows_cols(np.asarray([int(row)], dtype=np.int64), uniq)[0]
        t_pack1 = time.perf_counter()
        out = fetched[restore]
        self.stats.pack_seconds += time.perf_counter() - t_pack1
        return out
