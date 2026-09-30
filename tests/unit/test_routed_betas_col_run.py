"""RoutedBetas contiguous-column fast path must be byte-identical to the slow path.

A contiguous prefix request (the dense training idiom ``arr[:, :n]``) used to read
every pack column and discard most of them. When the needed pack columns form one
ascending contiguous run, a Zarr column slice touches only the covering column
chunks instead — measured ~6.5x faster on a 65 536-of-482 379 prefix. Correctness
is what this guards: same values, same NaN placement, including when a pack is
missing some columns (``col_map < 0``) so the fast path must *not* engage.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from mbs.matrix.virtual_hub_store import RoutedBetas


def _store(n_rows_a: int = 7, n_rows_b: int = 5, n_pack_cols: int = 40, n_loci: int = 24):
    """Two packs, interleaved virtual rows, so src_rows are scattered per pack."""
    rng = np.random.default_rng(0)
    a = rng.random((n_rows_a, n_pack_cols)).astype(np.float32)
    b = rng.random((n_rows_b, n_pack_cols)).astype(np.float32)
    a[0, 3] = np.nan  # keep a real NaN in play
    rows = []
    ia = ib = 0
    for v in range(n_rows_a + n_rows_b):
        if v % 2 == 0 and ia < n_rows_a:
            rows.append((v, "A", ia)); ia += 1
        elif ib < n_rows_b:
            rows.append((v, "B", ib)); ib += 1
        else:
            rows.append((v, "A", ia)); ia += 1
    route = pd.DataFrame(rows, columns=["row_index", "matrix_id", "src_row_index"])
    return a, b, route


def test_contiguous_prefix_matches_full_read() -> None:
    a, b, route = _store()
    n_loci = 24
    # identity col maps -> a prefix request maps to a contiguous pack-col run
    col_maps = {
        "A": np.arange(n_loci, dtype=np.int64),
        "B": np.arange(n_loci, dtype=np.int64),
    }
    rb = RoutedBetas(
        route=route, pack_arrays={"A": a, "B": b}, pack_col_maps=col_maps, n_loci=n_loci
    )
    got = rb[:, :12]
    # reference built independently from the pack arrays
    want = np.full((len(route), 12), np.nan, dtype=np.float32)
    for _, r in route.iterrows():
        src = a if r["matrix_id"] == "A" else b
        want[int(r["row_index"])] = src[int(r["src_row_index"]), :12]
    np.testing.assert_array_equal(np.isnan(got), np.isnan(want))
    np.testing.assert_allclose(got, want, equal_nan=True)


def test_missing_pack_columns_still_correct() -> None:
    """col_map < 0 (pack lacks a locus) must produce NaN, not a wrong slice."""
    a, b, route = _store()
    n_loci = 24
    col_maps = {
        "A": np.arange(n_loci, dtype=np.int64),
        "B": np.arange(n_loci, dtype=np.int64).copy(),
    }
    col_maps["B"][5] = -1  # pack B is missing virtual locus 5 -> breaks contiguity
    rb = RoutedBetas(
        route=route, pack_arrays={"A": a, "B": b}, pack_col_maps=col_maps, n_loci=n_loci
    )
    got = rb[:, :12]
    for _, r in route.iterrows():
        vrow = int(r["row_index"])
        if r["matrix_id"] == "B":
            assert np.isnan(got[vrow, 5]), "missing pack column must be NaN"
            keep = [c for c in range(12) if c != 5]
            np.testing.assert_allclose(
                got[vrow, keep], b[int(r["src_row_index"]), keep], equal_nan=True
            )
        else:
            np.testing.assert_allclose(
                got[vrow, :12], a[int(r["src_row_index"]), :12], equal_nan=True
            )


def test_non_contiguous_column_request_matches() -> None:
    """Fancy (non-contiguous) column selections must be unaffected by the fast path."""
    a, b, route = _store()
    n_loci = 24
    col_maps = {
        "A": np.arange(n_loci, dtype=np.int64),
        "B": np.arange(n_loci, dtype=np.int64),
    }
    rb = RoutedBetas(
        route=route, pack_arrays={"A": a, "B": b}, pack_col_maps=col_maps, n_loci=n_loci
    )
    cols = np.array([1, 7, 2, 19], dtype=np.int64)
    rows = np.array([0, 3, 4], dtype=np.int64)
    got = rb[np.ix_(rows, cols)]
    for i, vrow in enumerate(rows.tolist()):
        rec = route.loc[route["row_index"] == vrow].iloc[0]
        src = a if rec["matrix_id"] == "A" else b
        np.testing.assert_allclose(
            got[i], src[int(rec["src_row_index"]), cols], equal_nan=True
        )


def test_full_width_request_matches() -> None:
    """Requesting every column (fast path must not engage) still returns pack rows."""
    a, b, route = _store()
    n_loci = 24
    col_maps = {
        "A": np.arange(n_loci, dtype=np.int64),
        "B": np.arange(n_loci, dtype=np.int64),
    }
    rb = RoutedBetas(
        route=route, pack_arrays={"A": a, "B": b}, pack_col_maps=col_maps, n_loci=n_loci
    )
    got = rb[:, :n_loci]
    assert got.shape == (len(route), n_loci)
    for _, r in route.iterrows():
        src = a if r["matrix_id"] == "A" else b
        np.testing.assert_allclose(
            got[int(r["row_index"])], src[int(r["src_row_index"]), :n_loci], equal_nan=True
        )
