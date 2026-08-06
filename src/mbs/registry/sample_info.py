"""Export EWAS Data Hub sample-info archives to Parquet (prefer .txt over RData)."""

from __future__ import annotations

import zipfile
from pathlib import Path

import pandas as pd

# Primary phenotype column per family (Hub sample_*.txt).
FAMILY_VALUE_COLUMN: dict[str, str] = {
    "age": "age",
    "tissue": "tissue",
    "disease": "disease",
    "cancer": "disease",
    "blood": "cell_component",
    "brain": "tissue",
    "sex": "sex",
    "ancestry": "race",
    "bmi": "bmi",
}

# Backward-compatible private alias.
_FAMILY_VALUE_COLUMN = FAMILY_VALUE_COLUMN

_SAMPLE_ZIP_NAME: dict[str, str] = {
    "age": "sample_age_methylation_v1.zip",
    "tissue": "sample_tissue_methylation_v1.zip",
    "disease": "sample_disease_methylation_v1.zip",
    "cancer": "sample_cancer_methylation_v1.zip",
    "blood": "sample_blood_methylation_v1.zip",
    "brain": "sample_brain_methylation_v1.zip",
    "sex": "sample_sex_methylation_v1.zip",
    "ancestry": "sample_ancestry_category_methylation_v1.zip",
    "bmi": "sample_bmi_methylation_v1.zip",
}

_SAMPLE_TXT_STEM: dict[str, str] = {
    "age": "sample_age",
    "tissue": "sample_tissue",
    "disease": "sample_disease",
    "cancer": "sample_cancer",
    "blood": "sample_blood",
    "brain": "sample_brain",
    "sex": "sample_sex",
    "ancestry": "sample_ancestry_category",
    "bmi": "sample_bmi",
}


def sample_zip_filename(family: str) -> str:
    if family not in _SAMPLE_ZIP_NAME:
        raise ValueError(f"unsupported phenotype family: {family}")
    return _SAMPLE_ZIP_NAME[family]


def sample_txt_filename(family: str) -> str:
    if family not in _SAMPLE_TXT_STEM:
        raise ValueError(f"unsupported phenotype family: {family}")
    return f"{_SAMPLE_TXT_STEM[family]}.txt"


def unpacked_sample_info_dir(project_root: Path, family: str) -> Path:
    """Directory used by ``reports/inspection/ewas_datahub_samples/`` extracts."""
    if family not in _SAMPLE_ZIP_NAME:
        raise ValueError(f"unsupported phenotype family: {family}")
    zip_stem = Path(_SAMPLE_ZIP_NAME[family]).stem
    return project_root / "reports" / "inspection" / "ewas_datahub_samples" / zip_stem


def resolve_sample_info_txt(
    *,
    data_root: Path,
    family: str,
    project_root: Path | None = None,
    zip_path: Path | None = None,
    txt_path: Path | None = None,
) -> tuple[Path, str]:
    """Locate sample-info ``.txt`` from an explicit path, unpacked extract, or zip.

    Returns ``(path, source)`` where ``source`` is ``txt``, ``unpacked``, or ``zip``.
    """
    if family not in FAMILY_VALUE_COLUMN:
        raise ValueError(f"unsupported phenotype family: {family}")

    if txt_path is not None:
        resolved = txt_path.resolve()
        if not resolved.is_file():
            raise FileNotFoundError(f"sample-info txt not found: {resolved}")
        return resolved, "txt"

    if project_root is not None:
        unpacked = unpacked_sample_info_dir(project_root, family) / sample_txt_filename(family)
        if unpacked.is_file():
            return unpacked.resolve(), "unpacked"

    resolved_zip = (
        zip_path.resolve()
        if zip_path is not None
        else (data_root / "raw" / "ewas_datahub" / "download" / sample_zip_filename(family)).resolve()
    )
    if resolved_zip.is_file():
        return resolved_zip, "zip"

    raise FileNotFoundError(
        f"sample-info for family={family!r} not found as unpacked txt or zip "
        f"(looked under reports/inspection/ewas_datahub_samples/ and {resolved_zip})"
    )


def _find_txt_member(zf: zipfile.ZipFile) -> str:
    txt_members = [n for n in zf.namelist() if n.lower().endswith(".txt")]
    if not txt_members:
        raise FileNotFoundError(
            "sample-info zip has no .txt member; use scripts/export_ewas_sample_info.R "
            "for RData-only archives"
        )
    return sorted(txt_members)[0]


def read_r_style_table(path: Path) -> pd.DataFrame:
    """Read Hub sample_*.txt (R ``write.table`` with row names → use index_col=0)."""
    frame = pd.read_csv(
        path,
        sep=" ",
        quotechar='"',
        engine="python",
        na_values=["NA", "NaN"],
        dtype=str,
        index_col=0,
    )
    return frame.reset_index(drop=True)


# Backward-compatible private alias.
_read_r_style_table = read_r_style_table


def _normalize_sample_info_frame(frame: pd.DataFrame, *, family: str, source_name: str) -> pd.DataFrame:
    value_col = FAMILY_VALUE_COLUMN[family]
    if "sample_id" not in frame.columns:
        raise ValueError(f"sample-info missing sample_id column; got {list(frame.columns)[:12]}")
    if value_col not in frame.columns:
        raise ValueError(f"sample-info missing value column {value_col!r} for family={family}")

    keep = ["sample_id", "project_id", "platform", "sample_type", "sex", value_col]
    candidate_extra = ("tissue", "disease", "age", "bmi", "race", "cell_component")
    extra = [c for c in candidate_extra if c in frame.columns]
    cols: list[str] = []
    for c in keep + extra:
        if c in frame.columns and c not in cols:
            cols.append(c)
    out = frame.loc[:, cols].copy()
    rename_map = {
        "project_id": "study_id",
        value_col: "phenotype_value",
    }
    out.columns = [rename_map.get(str(c), str(c)) for c in out.columns]
    out["phenotype_family"] = family
    out["source_zip"] = source_name
    if family in {"age", "bmi"}:
        out["phenotype_value_numeric"] = pd.to_numeric(out["phenotype_value"], errors="coerce")
    return out


def export_sample_info_parquet(
    *,
    zip_path: Path | None = None,
    family: str,
    output_path: Path,
    txt_path: Path | None = None,
) -> Path:
    """Load sample-info from ``.txt`` or zip and write a normalized Parquet table."""
    if family not in FAMILY_VALUE_COLUMN:
        raise ValueError(f"unsupported phenotype family: {family}")

    output_path = output_path.resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if txt_path is not None:
        frame = read_r_style_table(txt_path.resolve())
        source_name = txt_path.name
    elif zip_path is not None:
        zip_path = zip_path.resolve()
        if not zip_path.is_file():
            raise FileNotFoundError(f"sample-info zip not found: {zip_path}")
        scratch = output_path.parent / f".extract_{family}"
        scratch.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(zip_path) as zf:
            member = _find_txt_member(zf)
            extracted = Path(zf.extract(member, path=scratch))
        frame = read_r_style_table(extracted)
        source_name = zip_path.name
    else:
        raise ValueError("export_sample_info_parquet requires zip_path or txt_path")

    out = _normalize_sample_info_frame(frame, family=family, source_name=source_name)
    out.to_parquet(output_path, index=False)
    return output_path


def export_family_from_data_root(
    data_root: Path,
    family: str,
    *,
    project_root: Path | None = None,
) -> Path:
    """Export one family using unpacked ``.txt`` when present, else download zip."""
    output = data_root / "canonical" / "phenotypes" / f"{family}_sample_info.parquet"
    resolved_root = project_root
    if resolved_root is None:
        # data_root is typically <project>/data
        candidate = data_root.resolve().parent
        if (candidate / "reports" / "inspection" / "ewas_datahub_samples").is_dir():
            resolved_root = candidate

    located, source = resolve_sample_info_txt(
        data_root=data_root,
        family=family,
        project_root=resolved_root,
    )
    if source == "zip":
        return export_sample_info_parquet(zip_path=located, family=family, output_path=output)
    return export_sample_info_parquet(txt_path=located, family=family, output_path=output)
