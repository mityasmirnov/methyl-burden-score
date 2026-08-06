"""Unified sample phenotype table + tissue ontology for Milestone 5c."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

PLATFORM_ALIASES = {
    "450K": "HM450",
    "450k": "HM450",
    "HM450": "HM450",
    "EPIC": "EPIC",
    "EPICv2": "EPICv2",
}

# Hub sex pack labels → binary class id (fail loud on anything else).
SEX_LABEL_TO_ID: dict[str, int] = {
    "Male": 0,
    "male": 0,
    "M": 0,
    "Female": 1,
    "female": 1,
    "F": 1,
}
SEX_CLASS_NAMES: tuple[str, ...] = ("Male", "Female")


def normalize_platform(raw: object | None) -> str | None:
    if raw is None or (isinstance(raw, float) and pd.isna(raw)):
        return None
    text = str(raw).strip()
    if text == "":
        return None
    return PLATFORM_ALIASES.get(text, text)


@dataclass(frozen=True, slots=True)
class TissueOntology:
    """Harmonized tissue class map (label → class_id)."""

    labels: tuple[str, ...]
    min_n: int
    label_to_id: dict[str, int]

    @property
    def class_names(self) -> list[str]:
        return list(self.labels)


def build_tissue_ontology(
    labels: list[str],
    *,
    min_n: int = 10,
    counts: dict[str, int] | None = None,
) -> TissueOntology:
    """Keep labels meeting min_n; assign stable sorted class ids."""
    if counts is None:
        counts = {}
        for lab in labels:
            counts[lab] = counts.get(lab, 0) + 1
    kept = sorted(lab for lab, n in counts.items() if n >= min_n)
    if not kept:
        raise ValueError(f"no tissue labels meet min_n={min_n}")
    return TissueOntology(
        labels=tuple(kept),
        min_n=min_n,
        label_to_id={lab: i for i, lab in enumerate(kept)},
    )


def tissue_ontology_to_dict(ontology: TissueOntology) -> dict[str, Any]:
    return {
        "version": "tissue-ontology-v1",
        "min_n": ontology.min_n,
        "classes": [
            {"class_id": ontology.label_to_id[lab], "label": lab} for lab in ontology.labels
        ],
    }


def write_tissue_ontology(path: Path, ontology: TissueOntology) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = yaml.safe_dump(tissue_ontology_to_dict(ontology), sort_keys=False)
    path.write_text(payload, encoding="utf-8")


def load_tissue_ontology(path: Path) -> TissueOntology:
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise TypeError(f"tissue ontology must be a mapping: {path}")
    classes = data.get("classes")
    if not isinstance(classes, list) or not classes:
        raise ValueError(f"tissue ontology missing classes: {path}")
    labels: list[str] = []
    label_to_id: dict[str, int] = {}
    for entry in classes:
        lab = str(entry["label"])
        cid = int(entry["class_id"])
        labels.append(lab)
        label_to_id[lab] = cid
    labels_sorted = tuple(sorted(labels, key=lambda x: label_to_id[x]))
    return TissueOntology(
        labels=labels_sorted,
        min_n=int(data.get("min_n", 10)),
        label_to_id=label_to_id,
    )


@dataclass(frozen=True, slots=True)
class SexOntology:
    """Binary sex class map (canonical Male/Female)."""

    labels: tuple[str, ...]
    label_to_id: dict[str, int]

    @property
    def class_names(self) -> list[str]:
        return list(self.labels)


def default_sex_ontology() -> SexOntology:
    return SexOntology(
        labels=SEX_CLASS_NAMES,
        label_to_id={"Male": 0, "Female": 1},
    )


def normalize_sex_label(raw: object | None) -> str | None:
    """Map Hub sex strings to canonical Male/Female; raise on unknown non-empty."""
    if raw is None or (isinstance(raw, float) and pd.isna(raw)):
        return None
    text = str(raw).strip()
    if text == "":
        return None
    if text not in SEX_LABEL_TO_ID:
        raise ValueError(f"unsupported sex label: {text!r}")
    return SEX_CLASS_NAMES[SEX_LABEL_TO_ID[text]]


def write_sex_ontology(path: Path, ontology: SexOntology | None = None) -> None:
    ont = ontology or default_sex_ontology()
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "version": "sex-ontology-v1",
        "classes": [
            {"class_id": ont.label_to_id[lab], "label": lab} for lab in ont.labels
        ],
    }
    path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")


def load_sex_ontology(path: Path) -> SexOntology:
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise TypeError(f"sex ontology must be a mapping: {path}")
    classes = data.get("classes")
    if not isinstance(classes, list) or not classes:
        raise ValueError(f"sex ontology missing classes: {path}")
    labels: list[str] = []
    label_to_id: dict[str, int] = {}
    for entry in classes:
        lab = str(entry["label"])
        cid = int(entry["class_id"])
        labels.append(lab)
        label_to_id[lab] = cid
    labels_sorted = tuple(sorted(labels, key=lambda x: label_to_id[x]))
    return SexOntology(labels=labels_sorted, label_to_id=label_to_id)


def _age_years(row: dict[str, Any]) -> float | None:
    for key in ("phenotype_value_numeric", "age"):
        raw = row.get(key)
        if raw is None or (isinstance(raw, float) and pd.isna(raw)):
            continue
        try:
            value = float(raw)
        except (TypeError, ValueError):
            continue
        if pd.isna(value):
            continue
        return value
    return None


def _tissue_label(row: dict[str, Any]) -> str | None:
    for key in ("phenotype_value", "tissue"):
        raw = row.get(key)
        if raw is None or (isinstance(raw, float) and pd.isna(raw)):
            continue
        text = str(raw).strip()
        if text == "":
            continue
        # Age pack stores numeric ages in phenotype_value; ignore those.
        try:
            float(text)
            continue
        except ValueError:
            return text
    return None


def _sex_from_row(row: dict[str, Any] | None, *, pack_primary: bool = False) -> str | None:
    """Extract canonical sex label from a sidecar/sample-info row.

    When ``pack_primary`` is True (sex-pack sidecar), try ``phenotype_value``
    then ``sex``. Otherwise only ``sex`` is consulted so numeric age
    ``phenotype_value`` is never treated as a sex label.
    """
    if row is None:
        return None
    keys = ("phenotype_value", "sex") if pack_primary else ("sex",)
    for key in keys:
        if key not in row:
            continue
        return normalize_sex_label(row.get(key))
    return None


def build_sample_phenotype_rows(
    *,
    age_pheno: pd.DataFrame,
    tissue_pheno: pd.DataFrame,
    sample_index: pd.DataFrame,
    matrix_id: str,
    ontology: TissueOntology,
    sex_pheno: pd.DataFrame | None = None,
    sex_ontology: SexOntology | None = None,
    source_system: str = "ewas_datahub",
) -> pd.DataFrame:
    """Build schema-shaped rows for samples present in ``sample_index``.

    ``age_pheno`` / ``tissue_pheno`` / optional ``sex_pheno`` are matrix sidecars
    or sample-info frames keyed by ``sample_id``. Samples need at least one of
    age / tissue / sex labels.
    """
    age_by = {str(r["sample_id"]): r for r in age_pheno.to_dict(orient="records")}
    tissue_by = {str(r["sample_id"]): r for r in tissue_pheno.to_dict(orient="records")}
    sex_by = (
        {str(r["sample_id"]): r for r in sex_pheno.to_dict(orient="records")}
        if sex_pheno is not None
        else {}
    )
    sex_ont = sex_ontology or default_sex_ontology()
    rows: list[dict[str, Any]] = []
    for rec in sample_index.to_dict(orient="records"):
        sid = str(rec["sample_id"])
        source_sid = str(rec.get("source_sample_id", sid))
        row_index = int(rec["row_index"])
        age_row = age_by.get(sid)
        tissue_row = tissue_by.get(sid)
        sex_row = sex_by.get(sid)

        age_years = _age_years(age_row) if age_row is not None else None
        if age_years is None and tissue_row is not None:
            age_years = _age_years(tissue_row)
        if age_years is None and sex_row is not None:
            age_years = _age_years(sex_row)
        age_mask = age_years is not None

        tissue_lab = _tissue_label(tissue_row) if tissue_row is not None else None
        if tissue_lab is None and age_row is not None:
            tissue_lab = _tissue_label(age_row)
        if tissue_lab is None and sex_row is not None:
            tissue_lab = _tissue_label(sex_row)
        tissue_class_id: int | None = None
        tissue_mask = False
        if tissue_lab is not None and tissue_lab in ontology.label_to_id:
            tissue_class_id = ontology.label_to_id[tissue_lab]
            tissue_mask = True
        elif tissue_lab is not None:
            tissue_lab = None  # filtered by min_n / ontology

        sex_lab = _sex_from_row(sex_row, pack_primary=True)
        if sex_lab is None:
            sex_lab = _sex_from_row(age_row)
        if sex_lab is None:
            sex_lab = _sex_from_row(tissue_row)
        sex_class_id: int | None = None
        sex_mask = False
        if sex_lab is not None:
            if sex_lab not in sex_ont.label_to_id:
                raise ValueError(f"sex label {sex_lab!r} not in sex ontology for {sid}")
            sex_class_id = sex_ont.label_to_id[sex_lab]
            sex_mask = True

        study_id = None
        platform_raw = None
        for src in (age_row, tissue_row, sex_row):
            if src is None:
                continue
            if study_id is None:
                study_id = src.get("study_id")
            if platform_raw is None:
                platform_raw = src.get("platform")
        if study_id is None or (isinstance(study_id, float) and pd.isna(study_id)):
            raise ValueError(f"missing study_id for sample_id={sid}")
        study_id = str(study_id)
        platform_id = normalize_platform(platform_raw)

        if not age_mask and not tissue_mask and not sex_mask:
            raise ValueError(f"sample {sid} has neither age, tissue, nor sex label")

        if age_mask and tissue_mask and sex_mask:
            family = "multi"
        elif sum(bool(x) for x in (age_mask, tissue_mask, sex_mask)) > 1:
            family = "multi"
        elif age_mask:
            family = "age"
        elif tissue_mask:
            family = "tissue"
        else:
            family = "sex"

        rows.append(
            {
                "sample_id": sid,
                "source_sample_id": source_sid,
                "study_id": study_id,
                "source_system": source_system,
                "phenotype_family": family,
                "platform_id": platform_id,
                "donor_id": study_id,
                "matrix_id": matrix_id,
                "row_index": row_index,
                "age_years": age_years,
                "age_mask": age_mask,
                "tissue_label": tissue_lab,
                "tissue_class_id": tissue_class_id,
                "tissue_mask": tissue_mask,
                "sex_label": sex_lab,
                "sex_class_id": sex_class_id,
                "sex_mask": sex_mask,
                "disease_mask": False,
                "cancer_mask": False,
                "blood_mask": False,
                "brain_mask": False,
                "sex": sex_lab,
            }
        )
    return pd.DataFrame(rows)


def write_sample_phenotype_table(path: Path, frame: pd.DataFrame) -> None:
    required = {
        "sample_id",
        "source_sample_id",
        "study_id",
        "source_system",
        "phenotype_family",
        "matrix_id",
        "row_index",
        "age_mask",
        "tissue_mask",
    }
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"sample phenotype table missing columns: {sorted(missing)}")
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(path, index=False)
