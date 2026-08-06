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


def build_sample_phenotype_rows(
    *,
    age_pheno: pd.DataFrame,
    tissue_pheno: pd.DataFrame,
    sample_index: pd.DataFrame,
    matrix_id: str,
    ontology: TissueOntology,
    source_system: str = "ewas_datahub",
) -> pd.DataFrame:
    """Build schema-shaped rows for samples present in ``sample_index``.

    ``age_pheno`` / ``tissue_pheno`` are matrix sidecars or sample-info frames
    keyed by ``sample_id``. Samples only in age get ``tissue_mask=False``;
    only in tissue get ``age_mask=False``; in both get both masks when labels
    resolve.
    """
    age_by = {str(r["sample_id"]): r for r in age_pheno.to_dict(orient="records")}
    tissue_by = {str(r["sample_id"]): r for r in tissue_pheno.to_dict(orient="records")}
    rows: list[dict[str, Any]] = []
    for rec in sample_index.to_dict(orient="records"):
        sid = str(rec["sample_id"])
        source_sid = str(rec.get("source_sample_id", sid))
        row_index = int(rec["row_index"])
        age_row = age_by.get(sid)
        tissue_row = tissue_by.get(sid)

        age_years = _age_years(age_row) if age_row is not None else None
        if age_years is None and tissue_row is not None:
            age_years = _age_years(tissue_row)
        age_mask = age_years is not None

        tissue_lab = _tissue_label(tissue_row) if tissue_row is not None else None
        if tissue_lab is None and age_row is not None:
            tissue_lab = _tissue_label(age_row)
        tissue_class_id: int | None = None
        tissue_mask = False
        if tissue_lab is not None and tissue_lab in ontology.label_to_id:
            tissue_class_id = ontology.label_to_id[tissue_lab]
            tissue_mask = True
        elif tissue_lab is not None:
            tissue_lab = None  # filtered by min_n / ontology

        study_id = None
        platform_raw = None
        sex = None
        if age_row is not None:
            study_id = age_row.get("study_id")
            platform_raw = age_row.get("platform")
            sex = age_row.get("sex")
        if tissue_row is not None:
            if study_id is None:
                study_id = tissue_row.get("study_id")
            if platform_raw is None:
                platform_raw = tissue_row.get("platform")
            if sex is None:
                sex = tissue_row.get("sex")
        if study_id is None or (isinstance(study_id, float) and pd.isna(study_id)):
            raise ValueError(f"missing study_id for sample_id={sid}")
        study_id = str(study_id)
        platform_id = normalize_platform(platform_raw)

        if age_mask and tissue_mask:
            family = "multi"
        elif age_mask:
            family = "age"
        elif tissue_mask:
            family = "tissue"
        else:
            raise ValueError(f"sample {sid} has neither age nor tissue label")

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
                "disease_mask": False,
                "cancer_mask": False,
                "blood_mask": False,
                "brain_mask": False,
                "sex": (
                    None if sex is None or (isinstance(sex, float) and pd.isna(sex)) else str(sex)
                ),
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
