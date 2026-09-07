"""Unit tests for RoutedBetas dense-style ``[:, :n]`` indexing (7H Phase 3)."""

from __future__ import annotations

from pathlib import Path
from uuid import uuid4

import numpy as np
import pandas as pd
import pytest

from mbs.matrix.store import (
    create_betas_zarr,
    matrix_store_paths,
    open_betas_zarr,
    read_sample_index,
    write_locus_index,
    write_sample_index,
)
from mbs.matrix.virtual_hub_store import (
    PACK_PRIORITY,
    VIRTUAL_MATRIX_ID,
    build_virtual_hub_store,
    open_betas_for_matrix,
)


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


@pytest.fixture
def workspace(monkeypatch: pytest.MonkeyPatch) -> Path:
    repo = _repo_root()
    scratch_base = repo / "scratch" / "pytest"
    scratch_base.mkdir(parents=True, exist_ok=True)
    root = scratch_base / f"7h-routed-{uuid4().hex}"
    root.mkdir()
    monkeypatch.setenv("MBS_ROOT", str(repo))
    monkeypatch.setenv("MBS_PROJECT_ROOT", str(repo))
    monkeypatch.setenv("MBS_DATA_ROOT", str(root / "data"))
    monkeypatch.setenv("MBS_SCRATCH_ROOT", str(root / "scratch"))
    monkeypatch.setenv("MBS_CACHE_ROOT", str(root / "cache"))
    monkeypatch.setenv("MBS_ARTIFACT_ROOT", str(root / "artifacts"))
    monkeypatch.setenv("MBS_DOCKER_ROOT", str(root / "docker"))
    return root


def _write_pack(
    data_root: Path,
    *,
    family: str,
    sample_ids: list[str],
    locus_ids: list[int],
    fill: float,
) -> Path:
    mid = f"matrix-hub-{family}-full-v1"
    root = data_root / "canonical" / "matrices" / mid
    root.mkdir(parents=True, exist_ok=True)
    write_sample_index(
        root / "sample_index.parquet",
        sample_ids=sample_ids,
        source_sample_ids=sample_ids,
    )
    n_loci = len(locus_ids)
    write_locus_index(
        root / "locus_index.parquet",
        locus_ids=np.asarray(locus_ids, dtype=np.uint64),
        canonical_keys=np.asarray([f"k{i}" for i in locus_ids], dtype=object),
        probe_ids=np.asarray([f"cg{i}" for i in locus_ids], dtype=object),
    )
    arr = create_betas_zarr(root / "betas.zarr", n_samples=len(sample_ids), n_loci=n_loci)
    for i in range(len(sample_ids)):
        arr[i, :] = np.full(n_loci, fill + 0.01 * i, dtype=np.float32)
    (root / "matrix_manifest.json").write_text(
        '{"matrix_id":"%s","shape":[%d,%d],"platform_id":"HM450"}'
        % (mid, len(sample_ids), n_loci),
        encoding="utf-8",
    )
    pheno = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "study_id": [f"GSE_{family}_{i}" for i in range(len(sample_ids))],
            "platform": ["450K"] * len(sample_ids),
            "phenotype_family": [family] * len(sample_ids),
        }
    )
    if family == "age":
        pheno["age"] = [30.0 + i for i in range(len(sample_ids))]
        pheno["phenotype_value_numeric"] = pheno["age"]
    elif family == "tissue":
        pheno["tissue"] = ["blood" if i % 2 == 0 else "brain" for i in range(len(sample_ids))]
    pheno.to_parquet(root / "sample_phenotypes.parquet", index=False)
    return root


def test_routed_betas_dense_column_slice_matches_row_vectors(workspace: Path) -> None:
    data_root = workspace / "data"
    _write_pack(
        data_root,
        family="age",
        sample_ids=["GSM_A", "GSM_BOTH"],
        locus_ids=[1, 2, 3],
        fill=0.1,
    )
    _write_pack(
        data_root,
        family="tissue",
        sample_ids=["GSM_T", "GSM_BOTH"],
        locus_ids=[2, 3, 4],
        fill=0.7,
    )
    result = build_virtual_hub_store(data_root=data_root)
    betas = open_betas_for_matrix(result.output_dir)
    assert betas.shape == (3, 2)

    dense = np.asarray(betas[:, :2], dtype=np.float32)
    assert dense.shape == (3, 2)
    for r in range(3):
        np.testing.assert_allclose(dense[r], betas[r], atol=1e-6)
        np.testing.assert_allclose(dense[r], betas[r, :], atol=1e-6)

    prefix = np.asarray(betas[:, :1], dtype=np.float32)
    assert prefix.shape == (3, 1)
    np.testing.assert_allclose(prefix[:, 0], dense[:, 0], atol=1e-6)

    sub = np.asarray(betas[1:3, :2], dtype=np.float32)
    np.testing.assert_allclose(sub, dense[1:3], atol=1e-6)


def test_routed_betas_matches_pack_source_and_index_invariant(workspace: Path) -> None:
    data_root = workspace / "data"
    _write_pack(
        data_root,
        family="age",
        sample_ids=["GSM_A", "GSM_BOTH"],
        locus_ids=[1, 2, 3],
        fill=0.1,
    )
    _write_pack(
        data_root,
        family="tissue",
        sample_ids=["GSM_T"],
        locus_ids=[2, 3, 4],
        fill=0.7,
    )
    result = build_virtual_hub_store(data_root=data_root)
    route = pd.read_parquet(result.route_path)
    betas = open_betas_for_matrix(result.output_dir)
    dense = np.asarray(betas[:, :], dtype=np.float32)

    sample_index = read_sample_index(matrix_store_paths(result.output_dir).sample_index_path)
    assert list(sample_index["row_index"].astype(int)) == list(range(len(sample_index)))
    assert list(route.sort_values("row_index")["row_index"].astype(int)) == list(
        range(len(route))
    )

    for rec in route.to_dict(orient="records"):
        row = int(rec["row_index"])
        mid = str(rec["matrix_id"])
        src_row = int(rec["src_row_index"])
        betas_path = Path(str(rec["betas_path"]))
        pack_arr = open_betas_zarr(
            betas_path if betas_path.name.endswith(".zarr") else betas_path / "betas.zarr"
        )
        assert pack_arr.shape[0] > src_row
        for col in range(betas.shape[1]):
            got = float(dense[row, col])
            via = float(betas[row, col])
            assert got == via or (np.isnan(got) and np.isnan(via))
        assert mid.startswith("matrix-hub-")
    assert PACK_PRIORITY[0] == "age"


def test_nine_pack_real_matrix_alignment_smoke() -> None:
    """Optional: ≥1 row from each pack × 64 cols vs pack Zarr (skip if absent)."""
    root = Path("/data/projects/methyl-burden-score/data/canonical/matrices") / VIRTUAL_MATRIX_ID
    if not (root / "route.parquet").is_file():
        pytest.skip("nine-pack virtual matrix not on disk")
    route = pd.read_parquet(root / "route.parquet")
    betas = open_betas_for_matrix(root)
    # One representative row per pack (all 9 families).
    rows: list[int] = []
    for mid, sub in route.groupby(route["matrix_id"].astype(str), sort=True):
        rows.append(int(sub["row_index"].iloc[0]))
    rows_a = np.asarray(rows, dtype=np.int64)
    cols = np.arange(64, dtype=np.int64)
    block = np.asarray(betas[np.ix_(rows_a, cols)], dtype=np.float32)
    via_slice = np.asarray(betas[rows_a, :64], dtype=np.float32)
    np.testing.assert_allclose(block, via_slice, equal_nan=True, atol=1e-5)
    assert len(rows_a) >= 8  # expect all nine packs when matrix is complete

    for i, r in enumerate(rows_a.tolist()):
        rec = route.loc[route["row_index"] == int(r)].iloc[0]
        mid = str(rec["matrix_id"])
        src_row = int(rec["src_row_index"])
        betas_path = Path(str(rec["betas_path"]))
        pack = open_betas_zarr(
            betas_path if betas_path.name.endswith(".zarr") else betas_path / "betas.zarr"
        )
        col_map = betas._pack_col_maps[mid]
        src_cols = col_map[cols]
        expected = np.full(64, np.nan, dtype=np.float32)
        valid = src_cols >= 0
        if valid.any():
            expected[valid] = np.asarray(
                pack[src_row, src_cols[valid].tolist()], dtype=np.float32
            )
        np.testing.assert_allclose(block[i], expected, equal_nan=True, atol=1e-5)
