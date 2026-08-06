"""Export EWAS Data Hub sample-info archives to Parquet (prefer .txt over RData)."""

from __future__ import annotations

import zipfile
from pathlib import Path

import pandas as pd

# Primary phenotype column per family (Hub sample_*.txt).
_FAMILY_VALUE_COLUMN: dict[str, str] = {
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


def sample_zip_filename(family: str) -> str:
    if family not in _SAMPLE_ZIP_NAME:
        raise ValueError(f"unsupported phenotype family: {family}")
    return _SAMPLE_ZIP_NAME[family]


def _find_txt_member(zf: zipfile.ZipFile) -> str:
    txt_members = [n for n in zf.namelist() if n.lower().endswith(".txt")]
    if not txt_members:
        raise FileNotFoundError(
            "sample-info zip has no .txt member; use scripts/export_ewas_sample_info.R "
            "for RData-only archives"
        )
    return sorted(txt_members)[0]


def _read_r_style_table(path: Path) -> pd.DataFrame:
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


def export_sample_info_parquet(
    *,
    zip_path: Path,
    family: str,
    output_path: Path,
) -> Path:
    """Extract sample-info .txt from zip and write a normalized Parquet table."""
    if family not in _FAMILY_VALUE_COLUMN:
        raise ValueError(f"unsupported phenotype family: {family}")
    zip_path = zip_path.resolve()
    if not zip_path.is_file():
        raise FileNotFoundError(f"sample-info zip not found: {zip_path}")

    scratch = output_path.parent / f".extract_{family}"
    scratch.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path) as zf:
        member = _find_txt_member(zf)
        extracted = Path(zf.extract(member, path=scratch))
    frame = _read_r_style_table(extracted)

    value_col = _FAMILY_VALUE_COLUMN[family]
    if "sample_id" not in frame.columns:
        raise ValueError(f"sample-info missing sample_id column; got {list(frame.columns)[:12]}")
    if value_col not in frame.columns:
        raise ValueError(f"sample-info missing value column {value_col!r} for family={family}")

    keep = ["sample_id", "project_id", "platform", "sample_type", "sex", value_col]
    # tissue pack uses tissue; age uses age; etc.
    candidate_extra = ("tissue", "disease", "age", "bmi", "race", "cell_component")
    extra = [c for c in candidate_extra if c in frame.columns]
    cols = []
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
    out["source_zip"] = zip_path.name
    # Coerce age/bmi numeric when possible
    if family in {"age", "bmi"}:
        out["phenotype_value_numeric"] = pd.to_numeric(out["phenotype_value"], errors="coerce")
    output_path = output_path.resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(output_path, index=False)
    return output_path


def export_family_from_data_root(data_root: Path, family: str) -> Path:
    zip_name = sample_zip_filename(family)
    zip_path = data_root / "raw" / "ewas_datahub" / "download" / zip_name
    output = data_root / "canonical" / "phenotypes" / f"{family}_sample_info.parquet"
    return export_sample_info_parquet(zip_path=zip_path, family=family, output_path=output)
