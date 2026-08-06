"""Phenotype labels for the GSE35069 pilot (CpGCorpus metadata join)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pyarrow as pa


@dataclass(frozen=True, slots=True)
class SamplePhenotype:
    sample_id: str
    cell_type: str
    donor_id: str
    title: str
    class_index: int


def _donor_from_title(title: str) -> str:
    """Parse donor suffix from titles like ``WB_1`` / ``CD4+_T_cells_3``."""
    if "_" not in title:
        raise ValueError(f"cannot parse donor from title without underscore: {title!r}")
    return title.rsplit("_", 1)[-1]


def load_gse35069_phenotypes(
    metadata_path: Path,
    *,
    sample_ids: list[str] | None = None,
) -> tuple[list[SamplePhenotype], list[str]]:
    """Load GSM→cell-type/donor labels from CpGCorpus ``metadata.arrow``.

    Returns ``(phenotypes ordered like sample_ids or arrow order, class_names)``.
    Raises if any requested sample_id is missing.
    """
    path = metadata_path.resolve()
    if not path.is_file():
        raise FileNotFoundError(f"phenotype metadata not found: {path}")

    table = pa.ipc.open_file(path).read_all()
    required = {"GSM_ID", "title", "tissue/cell type:ch1"}
    missing_cols = required - set(table.column_names)
    if missing_cols:
        raise ValueError(f"metadata.arrow missing columns: {sorted(missing_cols)}")

    gsms = [str(x) for x in table.column("GSM_ID").to_pylist()]
    titles = [str(x) for x in table.column("title").to_pylist()]
    cell_types = [str(x) for x in table.column("tissue/cell type:ch1").to_pylist()]

    by_gsm: dict[str, tuple[str, str]] = {}
    for gsm, title, cell in zip(gsms, titles, cell_types, strict=True):
        by_gsm[gsm] = (cell, title)

    ordered_ids = sample_ids if sample_ids is not None else gsms
    class_names = sorted({by_gsm[sid][0] for sid in ordered_ids if sid in by_gsm})
    if not class_names:
        raise ValueError("no overlapping sample IDs between matrix and metadata")
    class_to_idx = {name: i for i, name in enumerate(class_names)}

    phenotypes: list[SamplePhenotype] = []
    missing: list[str] = []
    for sid in ordered_ids:
        if sid not in by_gsm:
            missing.append(sid)
            continue
        cell, title = by_gsm[sid]
        phenotypes.append(
            SamplePhenotype(
                sample_id=sid,
                cell_type=cell,
                donor_id=_donor_from_title(title),
                title=title,
                class_index=class_to_idx[cell],
            )
        )
    if missing:
        raise KeyError(f"metadata missing {len(missing)} sample_id(s); first={missing[0]!r}")
    return phenotypes, class_names
