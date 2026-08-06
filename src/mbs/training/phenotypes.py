"""Phenotype labels for pilot and Hub sample-info joins."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import pyarrow as pa


@dataclass(frozen=True, slots=True)
class SamplePhenotype:
    sample_id: str
    cell_type: str
    donor_id: str
    title: str
    class_index: int
    study_id: str | None = None
    age: float | None = None
    platform: str | None = None


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
                study_id="GSE35069",
            )
        )
    if missing:
        raise KeyError(f"metadata missing {len(missing)} sample_id(s); first={missing[0]!r}")
    return phenotypes, class_names


def load_hub_sample_info_phenotypes(
    parquet_path: Path,
    *,
    sample_ids: list[str] | None = None,
    value_column: str = "phenotype_value",
    empty_as_control: bool = False,
) -> tuple[list[SamplePhenotype], list[str]]:
    """Load multiclass labels from Hub sample-info / matrix sidecar Parquet.

    When ``empty_as_control`` is True, blank/NA ``value_column`` values become
    the class ``control`` (disease case/control packs).
    """
    path = parquet_path.resolve()
    if not path.is_file():
        raise FileNotFoundError(f"sample-info parquet not found: {path}")
    frame = pd.read_parquet(path)
    if "sample_id" not in frame.columns:
        raise ValueError("parquet must contain sample_id")
    if value_column not in frame.columns:
        raise ValueError(f"parquet must contain {value_column}")
    records = frame.to_dict(orient="records")
    by_id: dict[str, dict[str, object]] = {str(row["sample_id"]): row for row in records}
    ordered = sample_ids if sample_ids is not None else [str(x) for x in frame["sample_id"]]

    def _label_of(row: dict[str, object]) -> str | None:
        raw = row.get(value_column)
        if raw is None or (isinstance(raw, float) and pd.isna(raw)):
            return "control" if empty_as_control else None
        text = str(raw).strip()
        if text == "":
            return "control" if empty_as_control else None
        return text

    labels = []
    for sid in ordered:
        if sid not in by_id:
            continue
        lab = _label_of(by_id[sid])
        if lab is not None:
            labels.append(lab)
    class_names = sorted(set(labels))
    if not class_names:
        raise ValueError("no phenotype labels found in sample-info parquet")
    class_to_idx = {name: i for i, name in enumerate(class_names)}
    phenotypes: list[SamplePhenotype] = []
    for sid in ordered:
        if sid not in by_id:
            raise KeyError(f"sample_id missing from sample-info: {sid}")
        row = by_id[sid]
        value = _label_of(row)
        if value is None:
            raise KeyError(f"missing phenotype for sample_id={sid}")
        study = row.get("study_id")
        age_val = row.get("phenotype_value_numeric")
        if age_val is None and "age" in row:
            age_val = row.get("age")
        if age_val is not None:
            try:
                age_f = float(age_val)  # type: ignore[arg-type]
                if pd.isna(age_f):
                    age_f = None
            except (TypeError, ValueError):
                age_f = None
        else:
            age_f = None
        platform = row.get("platform")
        phenotypes.append(
            SamplePhenotype(
                sample_id=sid,
                cell_type=str(value),
                donor_id=str(study or sid),
                title=str(value),
                class_index=class_to_idx[str(value)],
                study_id=None if study is None else str(study),
                age=age_f,
                platform=None if platform is None else str(platform),
            )
        )
    return phenotypes, class_names


def load_hub_regression_phenotypes(
    parquet_path: Path,
    *,
    sample_ids: list[str] | None = None,
    value_column: str = "phenotype_value_numeric",
) -> tuple[list[SamplePhenotype], list[str]]:
    """Load continuous targets (age/BMI) from Hub sidecar / sample-info Parquet."""
    path = parquet_path.resolve()
    if not path.is_file():
        raise FileNotFoundError(f"sample-info parquet not found: {path}")
    frame = pd.read_parquet(path)
    if "sample_id" not in frame.columns:
        raise ValueError("parquet must contain sample_id")
    if value_column not in frame.columns and value_column == "phenotype_value_numeric":
        if "age" in frame.columns:
            value_column = "age"
        else:
            raise ValueError("parquet must contain phenotype_value_numeric or age")
    records = frame.to_dict(orient="records")
    by_id: dict[str, dict[str, object]] = {str(row["sample_id"]): row for row in records}
    ordered = sample_ids if sample_ids is not None else [str(x) for x in frame["sample_id"]]
    phenotypes: list[SamplePhenotype] = []
    for sid in ordered:
        if sid not in by_id:
            raise KeyError(f"sample_id missing from sample-info: {sid}")
        row = by_id[sid]
        raw = row.get(value_column)
        try:
            age_f = float(raw)  # type: ignore[arg-type]
        except (TypeError, ValueError):
            raise KeyError(f"non-numeric regression target for sample_id={sid}") from None
        if pd.isna(age_f):
            raise KeyError(f"missing regression target for sample_id={sid}")
        study = row.get("study_id")
        platform = row.get("platform")
        phenotypes.append(
            SamplePhenotype(
                sample_id=sid,
                cell_type="age",
                donor_id=str(study or sid),
                title=str(age_f),
                class_index=0,
                study_id=None if study is None else str(study),
                age=age_f,
                platform=None if platform is None else str(platform),
            )
        )
    return phenotypes, ["age"]
