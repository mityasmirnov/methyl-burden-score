"""Versioned deepmat-data-v1 release: catalog population, census, eligibility."""

from __future__ import annotations

import hashlib
import json
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast
from urllib.error import URLError
from urllib.request import urlopen

import duckdb
import pandas as pd

from mbs.annotation.build import DEFAULT_GRAPH_ID
from mbs.annotation.manifest import git_commit, sha256_file, utc_now_iso, write_json
from mbs.atlas_study_enrichment import (
    build_study_atlas_enrichment,
    merge_atlas_enrichment_into_studies,
    write_study_atlas_enrichment_report,
)
from mbs.catalog import build_catalog
from mbs.geo_metadata import (
    GEO_PARQUET_NAME,
    geo_backfill_enabled,
    load_census_snapshot,
    load_geo_frame,
    merge_geo_sample_metadata,
    write_geo_backfill_pilot_report,
)
from mbs.paths import DataPaths
from mbs.platform_id import PLATFORM_ALIASES, normalize_platform
from mbs.registry.sample_info import FAMILY_VALUE_COLUMN
from mbs.static_features.export_cpgpt import DEFAULT_FEATURE_SET_ID

RELEASE_ID = "deepmat-data-v1"
ARTIFACT_VERSION = "release-manifest-v1"
EWAS_DB_ADVERTISED_N = 1989
EWAS_DB_INDEX_URL = "https://download.cncb.ac.cn/ewas/datahub/EWAS_db/"
HUB_SOURCE_RELEASE_ID = "ewas-datahub-baseline-v1"
EWAS_DB_SOURCE_RELEASE_ID = "ewas-datahub-db-v1"
FROZEN_SPLIT_RUN_ID = "stage0-flat-deeprvat-age-tissue-sex-full-v1"
_COMMIT_RE = re.compile(r"^[0-9a-f]{40}$")
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")

HUB_FAMILIES: tuple[str, ...] = tuple(FAMILY_VALUE_COLUMN.keys())

# Programme cutoffs (docs/plans/post-v0-scientific-programme.md).
AGE_BMI_MIN_SAMPLES = 1000
AGE_BMI_MIN_STUDIES = 5
DISEASE_MIN_CASES = 200
DISEASE_MIN_CONTROLS = 200
DISEASE_MIN_STUDIES = 3
TISSUE_MIN_PER_CLASS = 100
TISSUE_MIN_STUDIES_PER_CLASS = 2
SEX_MIN_SAMPLES = 200

_CONTROL_SAMPLE_TYPES = frozenset(
    {
        "control",
        "controls",
        "healthy",
        "normal",
        "unaffected",
    }
)


def _as_float(value: object) -> float | None:
    """Coerce a scalar to float; None if missing/non-numeric."""
    if _is_blank(value):
        return None
    number = pd.to_numeric(value, errors="coerce")  # type: ignore[arg-type]
    if number is None or (isinstance(number, float) and pd.isna(number)):
        return None
    try:
        return float(number)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None


def _as_int(value: object) -> int | None:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    return int(value)  # type: ignore[arg-type]


def _as_dataframe(frame: pd.DataFrame | pd.Series) -> pd.DataFrame:
    if isinstance(frame, pd.Series):
        return frame.to_frame().T
    return frame


@dataclass(frozen=True, slots=True)
class ReleasePaths:
    """Filesystem layout for a versioned data release."""

    root: Path
    manifest_path: Path
    catalog_db: Path
    catalog_tables: Path
    phenotypes_dir: Path
    matrices_dir: Path
    ontologies_dir: Path
    splits_dir: Path


@dataclass(frozen=True, slots=True)
class RefreshResult:
    release_id: str
    release_root: str
    catalog_path: str
    n_samples: int
    n_studies: int
    n_phenotype_rows: int
    ewas_db_n_local_studies: int
    ewas_db_n_local_gsm: int
    ewas_db_mirror_complete: bool
    report_dir: str | None


def release_paths(data_root: Path, release_id: str = RELEASE_ID) -> ReleasePaths:
    root = (data_root / "canonical" / "releases" / release_id).resolve()
    catalog = root / "catalog"
    return ReleasePaths(
        root=root,
        manifest_path=root / "release_manifest.json",
        catalog_db=catalog / "catalog.duckdb",
        catalog_tables=catalog / "tables",
        phenotypes_dir=root / "phenotypes",
        matrices_dir=root / "matrices",
        ontologies_dir=root / "ontologies",
        splits_dir=root / "splits",
    )


def release_inspection_report_dir(project_root: Path, release_id: str = RELEASE_ID) -> Path:
    """Inspection reports use underscore slug (``deepmat_data_v1``), not release id hyphen."""
    return project_root / "reports" / "inspection" / release_id.replace("-", "_")


def path_size_digest(relative_path: str, byte_size: int) -> str:
    """Content-address of path+size (no full-file hash)."""
    return hashlib.sha256(f"{relative_path}:{byte_size}".encode()).hexdigest()


def _is_blank(value: object) -> bool:
    if value is None:
        return True
    if isinstance(value, float) and pd.isna(value):
        return True
    text = str(value).strip()
    return text == "" or text.lower() in {"nan", "none", "na", "<na>"}


def _is_control_sample_type(value: object) -> bool:
    if _is_blank(value):
        return False
    text = str(value).strip().lower()
    if text in _CONTROL_SAMPLE_TYPES:
        return True
    return text.startswith("control")


def scan_ewas_db_tree(ewas_db_root: Path) -> pd.DataFrame:
    """Shallow inventory of ``EWAS_db/{STUDY}/GSM*.txt`` (no beta reads)."""
    ewas_db_root = ewas_db_root.resolve()
    rows: list[dict[str, Any]] = []
    if not ewas_db_root.is_dir():
        return pd.DataFrame(
            columns=[
                "study_id",
                "sample_id",
                "path",
                "byte_size",
                "mtime",
                "sha256",
            ]
        )
    for study_dir in sorted(p for p in ewas_db_root.iterdir() if p.is_dir()):
        study_id = study_dir.name
        for path in sorted(p for p in study_dir.iterdir() if p.is_file() and p.suffix == ".txt"):
            try:
                stat = path.stat()
            except OSError:
                continue
            rel = path.relative_to(ewas_db_root).as_posix()
            size = int(stat.st_size)
            rows.append(
                {
                    "study_id": study_id,
                    "sample_id": path.stem,
                    "path": str(path.resolve()),
                    "byte_size": size,
                    "mtime": float(stat.st_mtime),
                    "sha256": path_size_digest(rel, size),
                }
            )
    return pd.DataFrame(rows)


def ewas_db_study_inventory(assay_files: pd.DataFrame) -> pd.DataFrame:
    if assay_files.empty:
        return pd.DataFrame(columns=["study_id", "n_txt", "total_bytes", "newest_mtime"])
    grouped = assay_files.groupby("study_id", as_index=False).agg(
        n_txt=("sample_id", "count"),
        total_bytes=("byte_size", "sum"),
        newest_mtime=("mtime", "max"),
    )
    return cast(pd.DataFrame, grouped).sort_values(by=["study_id"])


def fetch_remote_ewas_db_study_names(*, timeout_s: float = 120.0) -> list[str] | None:
    """One HTTP GET of the EWAS_db index; return study names or None on failure."""
    try:
        with urlopen(EWAS_DB_INDEX_URL, timeout=timeout_s) as response:
            html = response.read().decode("utf-8", errors="replace")
    except (URLError, TimeoutError, OSError):
        return None
    names = sorted(
        {
            match.rstrip("/")
            for match in re.findall(r'href="([^"]+/)"', html)
            if match not in {"../", "/"}
            and not match.startswith(("http", "javascript", "#"))
            and "(" not in match
        }
    )
    return [n for n in names if n and n != "index.html"]


def _label_status_for_row(
    *,
    family: str,
    phenotype_id: str,
    value: object,
    sample_type: object,
) -> tuple[str, bool]:
    """Return (label_status, is_observed)."""
    if _is_blank(value):
        return "unknown", False
    if phenotype_id in {"disease", "cancer"}:
        if _is_control_sample_type(sample_type):
            return "control", True
        return "case", True
    if family in {"disease", "cancer"} and phenotype_id == family:
        if _is_control_sample_type(sample_type):
            return "control", True
        return "case", True
    return "observed", True


def _donor_replicate_from_row(row: dict[str, Any]) -> tuple[str | None, str | None]:
    """Copy donor/replicate when Hub sample-info provides them; never invent."""
    donor: str | None = None
    for key in ("donor_id", "subject_id", "individual_id"):
        raw = row.get(key)
        if raw is None or (isinstance(raw, float) and pd.isna(raw)):
            continue
        text = str(raw).strip()
        if text:
            donor = text
            break
    replicate: str | None = None
    for key in ("replicate_group", "replicate_id", "technical_replicate"):
        raw = row.get(key)
        if raw is None or (isinstance(raw, float) and pd.isna(raw)):
            continue
        text = str(raw).strip()
        if text:
            replicate = text
            break
    return donor, replicate


def build_long_form_phenotypes(data_root: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build long-form phenotypes, membership, and sample stubs from Hub packs."""
    pheno_rows: list[dict[str, Any]] = []
    membership_rows: list[dict[str, Any]] = []
    sample_meta: dict[str, dict[str, Any]] = {}

    for family in HUB_FAMILIES:
        path = data_root / "canonical" / "phenotypes" / f"{family}_sample_info.parquet"
        if not path.is_file():
            continue
        frame = pd.read_parquet(path)
        if "sample_id" not in frame.columns or "study_id" not in frame.columns:
            raise ValueError(f"{path} missing sample_id/study_id")
        primary_pid = family if family != "ancestry" else "ancestry"
        for row_index, row in enumerate(frame.to_dict(orient="records")):
            sid = str(row["sample_id"]).strip()
            study_id = str(row["study_id"]).strip()
            if not sid or not study_id or _is_blank(sid) or _is_blank(study_id):
                raise ValueError(f"blank sample/study in {family} row {row_index}")
            platform_id = normalize_platform(row.get("platform"))
            sample_type = row.get("sample_type")
            donor_id, replicate_group = _donor_replicate_from_row(row)
            meta = sample_meta.setdefault(
                sid,
                {
                    "sample_id": sid,
                    "source_sample_id": sid,
                    "study_id": study_id,
                    "platform_id": platform_id,
                    "age": None,
                    "sex": None,
                    "tissue_raw": None,
                    "case_control": None,
                    "donor_id": donor_id,
                    "replicate_group": replicate_group,
                },
            )
            if meta["platform_id"] is None and platform_id is not None:
                meta["platform_id"] = platform_id
            if meta.get("donor_id") is None and donor_id is not None:
                meta["donor_id"] = donor_id
            if meta.get("replicate_group") is None and replicate_group is not None:
                meta["replicate_group"] = replicate_group

            membership_rows.append(
                {
                    "sample_id": sid,
                    "source_release_id": HUB_SOURCE_RELEASE_ID,
                    "phenotype_family": family,
                    "source_file": path.name,
                    "source_row": row_index,
                    "matrix_id": None,
                    "row_index": None,
                }
            )

            observations: list[tuple[str, object]] = []
            if "phenotype_value" in row:
                observations.append((primary_pid, row.get("phenotype_value")))
            if "phenotype_value_numeric" in row and primary_pid in {"age", "bmi"}:
                observations = [(primary_pid, row.get("phenotype_value_numeric"))]
            for col, pid in (
                ("sex", "sex"),
                ("age", "age"),
                ("tissue", "tissue"),
                ("bmi", "bmi"),
            ):
                if col in row and pid != primary_pid:
                    observations.append((pid, row.get(col)))

            seen_pids: set[str] = set()
            for phenotype_id, value in observations:
                if phenotype_id in seen_pids:
                    continue
                seen_pids.add(phenotype_id)
                label_status, is_observed = _label_status_for_row(
                    family=family,
                    phenotype_id=phenotype_id,
                    value=value,
                    sample_type=sample_type,
                )
                numeric_value = None
                categorical_value = None
                if phenotype_id in {"age", "bmi"} and not _is_blank(value):
                    numeric_value = _as_float(value)
                    if numeric_value is None:
                        label_status, is_observed = "unknown", False
                elif not _is_blank(value):
                    categorical_value = str(value).strip()

                pheno_rows.append(
                    {
                        "sample_id": sid,
                        "phenotype_id": phenotype_id,
                        "numeric_value": numeric_value,
                        "categorical_value": categorical_value,
                        "label_status": label_status,
                        "is_observed": is_observed,
                        "source_family": family,
                        "source_record_id": f"{family}:{row_index}",
                        "ontology_id": None,
                    }
                )

                if is_observed and phenotype_id == "age" and numeric_value is not None:
                    meta["age"] = numeric_value
                if is_observed and phenotype_id == "sex" and categorical_value is not None:
                    meta["sex"] = categorical_value
                if is_observed and phenotype_id == "tissue" and categorical_value is not None:
                    meta["tissue_raw"] = categorical_value
                if phenotype_id in {"disease", "cancer"} and label_status in {
                    "case",
                    "control",
                }:
                    meta["case_control"] = label_status

    phenotypes = pd.DataFrame(pheno_rows)
    membership = pd.DataFrame(membership_rows)
    samples = pd.DataFrame(list(sample_meta.values())) if sample_meta else pd.DataFrame()
    if not phenotypes.empty:
        # Hub packs can repeat GSM within a family; long-form PK is unique per family.
        phenotypes = phenotypes.drop_duplicates(
            subset=["sample_id", "phenotype_id", "source_family"],
            keep="first",
        )
    if not membership.empty:
        membership = membership.drop_duplicates(
            subset=["sample_id", "source_release_id", "phenotype_family"],
            keep="first",
        )
    return phenotypes, membership, samples


def _seed_phenotype_table() -> pd.DataFrame:
    rows = [
        ("age", "Age", "continuous", None, "years"),
        ("bmi", "BMI", "continuous", None, "kg/m2"),
        ("tissue", "Tissue", "categorical", None, None),
        ("sex", "Sex", "categorical", None, None),
        ("disease", "Disease", "categorical", None, None),
        ("cancer", "Cancer", "categorical", None, None),
        ("blood", "Blood cell type", "categorical", None, None),
        ("brain", "Brain region", "categorical", None, None),
        ("ancestry", "Ancestry category", "categorical", None, None),
    ]
    return pd.DataFrame(
        rows,
        columns=["phenotype_id", "phenotype_name", "phenotype_type", "ontology_id", "unit"],
    )


def _seed_platforms() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "platform_id": "HM450",
                "platform_name": "Illumina HumanMethylation450",
                "manufacturer": "Illumina",
                "manifest_version": None,
                "nominal_probe_count": 485577,
                "genome_build": "GRCh38",
            },
            {
                "platform_id": "EPIC",
                "platform_name": "Illumina MethylationEPIC",
                "manufacturer": "Illumina",
                "manifest_version": None,
                "nominal_probe_count": 866836,
                "genome_build": "GRCh38",
            },
            {
                "platform_id": "EPICv2",
                "platform_name": "Illumina MethylationEPICv2",
                "manufacturer": "Illumina",
                "manifest_version": None,
                "nominal_probe_count": None,
                "genome_build": "GRCh38",
            },
        ]
    )


def scan_matrix_artifacts(data_root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Index canonical matrix manifests (no Zarr reads)."""
    matrices_root = data_root / "canonical" / "matrices"
    artifact_rows: list[dict[str, Any]] = []
    sample_rows: list[dict[str, Any]] = []
    if not matrices_root.is_dir():
        return pd.DataFrame(), pd.DataFrame()
    for matrix_dir in sorted(p for p in matrices_root.iterdir() if p.is_dir()):
        manifest_path = matrix_dir / "matrix_manifest.json"
        if not manifest_path.is_file():
            continue
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        matrix_id = str(manifest.get("matrix_id") or matrix_dir.name)
        shape = manifest.get("shape") or [None, None]
        family = None
        if isinstance(manifest.get("notes"), str) and "pack" in manifest["notes"]:
            family = None
        notes = manifest.get("notes")
        if isinstance(notes, str):
            for fam in HUB_FAMILIES:
                hit = (
                    f" {fam};" in f" {notes}"
                    or f"pack {fam}" in notes
                    or f"baseline pack {fam}" in notes
                )
                if hit:
                    family = fam
                    break
        n_loci = None
        n_samples_m: int | None = None
        if isinstance(shape, list) and len(shape) >= 2:
            n_samples_m = _as_int(shape[0])
            n_loci = _as_int(shape[1])
        artifact_rows.append(
            {
                "matrix_id": matrix_id,
                "path": str(matrix_dir.resolve()),
                "platform_id": manifest.get("platform_id"),
                "processing_level": manifest.get("processing_level"),
                "genome_build": manifest.get("genome_build"),
                "n_samples": n_samples_m,
                "n_loci": n_loci,
                "manifest_sha256": sha256_file(manifest_path),
                "phenotype_family": family,
                "notes": notes,
            }
        )
        sample_index = matrix_dir / "sample_index.parquet"
        if sample_index.is_file():
            idx = pd.read_parquet(sample_index)
            if "sample_id" not in idx.columns:
                raise ValueError(f"sample_index missing sample_id: {sample_index}")
            for i, rec in enumerate(idx.to_dict(orient="records")):
                sample_rows.append(
                    {
                        "matrix_id": matrix_id,
                        "sample_id": str(rec["sample_id"]),
                        "row_index": int(rec.get("row_index", i)),
                        "source_sample_id": str(rec.get("source_sample_id") or rec["sample_id"]),
                    }
                )
    return pd.DataFrame(artifact_rows), pd.DataFrame(sample_rows)


def load_frozen_split(artifact_root: Path) -> pd.DataFrame:
    """Ingest frozen 5d study-grouped split into fold_assignment rows."""
    split_path = artifact_root / "runs" / FROZEN_SPLIT_RUN_ID / "split.json"
    if not split_path.is_file():
        return pd.DataFrame(
            columns=["split_id", "sample_id", "outer_fold", "inner_fold", "group_id", "role"]
        )
    payload = json.loads(split_path.read_text(encoding="utf-8"))
    split_id = str(payload.get("split_id") or FROZEN_SPLIT_RUN_ID)
    rows = [
        {
            "split_id": split_id,
            "sample_id": str(sample["sample_id"]),
            "outer_fold": 0,
            "inner_fold": None,
            "group_id": str(sample.get("group_id") or sample.get("study_id")),
            "role": str(sample["role"]),
        }
        for sample in payload.get("samples") or []
    ]
    return pd.DataFrame(rows)


def compute_trait_eligibility(
    phenotypes: pd.DataFrame,
    samples: pd.DataFrame,
) -> pd.DataFrame:
    """Apply programme cutoffs; unknown disease labels are never controls."""
    if phenotypes.empty:
        return pd.DataFrame(
            columns=[
                "phenotype_id",
                "phenotype_family",
                "task_type",
                "n_samples",
                "n_cases",
                "n_controls",
                "n_unknown",
                "prevalence",
                "n_studies",
                "n_platforms",
                "n_tissues",
                "eligible_core_task",
                "eligible_auxiliary_task",
                "eligible_external_evaluation",
                "exclusion_reason",
            ]
        )
    study_by_sample = samples.set_index("sample_id")["study_id"].to_dict()
    if "platform_id" in samples.columns:
        platform_by_sample = samples.set_index("sample_id")["platform_id"].to_dict()
    else:
        platform_by_sample = {}
    if "tissue_raw" in samples.columns:
        tissue_by_sample = samples.set_index("sample_id")["tissue_raw"].to_dict()
    else:
        tissue_by_sample = {}

    rows: list[dict[str, Any]] = []
    grouped = phenotypes.groupby(["phenotype_id", "source_family"], sort=True)
    for key, group in grouped:
        phenotype_id, source_family = cast(tuple[str, str], key)
        observed = _as_dataframe(group[group["is_observed"]])
        sample_ids = set(observed["sample_id"].astype(str))
        n_samples = len(sample_ids)
        studies = {study_by_sample.get(sid) for sid in sample_ids}
        studies.discard(None)
        platforms = {platform_by_sample.get(sid) for sid in sample_ids}
        platforms.discard(None)
        n_studies = len(studies)
        n_platforms = len(platforms)

        n_cases = int((observed["label_status"] == "case").sum())
        n_controls = int((observed["label_status"] == "control").sum())
        n_unknown = int((group["label_status"] == "unknown").sum())
        # Unique samples for case/control when multi-rows
        if phenotype_id in {"disease", "cancer"}:
            n_cases = int(
                pd.Series(observed.loc[observed["label_status"] == "case", "sample_id"]).nunique()
            )
            n_controls = int(
                pd.Series(
                    observed.loc[observed["label_status"] == "control", "sample_id"]
                ).nunique()
            )
            n_unknown = int(
                pd.Series(group.loc[group["label_status"] == "unknown", "sample_id"]).nunique()
            )

        prevalence = None
        if phenotype_id in {"disease", "cancer"} and (n_cases + n_controls) > 0:
            prevalence = float(n_cases) / float(n_cases + n_controls)

        tissues = {
            tissue_by_sample[sid]
            for sid in sample_ids
            if sid in tissue_by_sample and not _is_blank(tissue_by_sample[sid])
        }

        task_type = "categorical"
        if phenotype_id in {"age", "bmi"}:
            task_type = "continuous"
        elif phenotype_id in {"disease", "cancer"}:
            task_type = "binary_or_multilabel"
        elif phenotype_id == "sex":
            task_type = "binary"
        elif phenotype_id in {"tissue", "blood", "brain", "ancestry"}:
            task_type = "multiclass"

        eligible_core = False
        eligible_aux = False
        eligible_ext = False
        reason: str | None = None

        if phenotype_id in {"age", "bmi"}:
            ages = pd.Series(observed["numeric_value"]).dropna()
            multi_study_range = False
            if n_studies > 1 and not ages.empty:
                multi_study_range = True
            eligible_core = (
                n_samples >= AGE_BMI_MIN_SAMPLES
                and n_studies >= AGE_BMI_MIN_STUDIES
                and multi_study_range
            )
            if not eligible_core:
                reason = (
                    f"need ≥{AGE_BMI_MIN_SAMPLES} samples, ≥{AGE_BMI_MIN_STUDIES} studies, "
                    "range across >1 study"
                )
            eligible_aux = n_samples >= 100 and n_studies >= 2
        elif phenotype_id in {"disease", "cancer"}:
            eligible_core = (
                n_cases >= DISEASE_MIN_CASES
                and n_controls >= DISEASE_MIN_CONTROLS
                and n_studies >= DISEASE_MIN_STUDIES
            )
            if not eligible_core:
                reason = (
                    f"need ≥{DISEASE_MIN_CASES} cases, ≥{DISEASE_MIN_CONTROLS} controls, "
                    f"≥{DISEASE_MIN_STUDIES} studies (unknown≠control)"
                )
            eligible_aux = n_cases >= 100 and n_studies >= 2
            eligible_ext = n_cases >= 50
        elif phenotype_id == "tissue":
            class_counts = (
                observed.groupby("categorical_value")["sample_id"].nunique()
                if not observed.empty
                else pd.Series(dtype=int)
            )
            strong = class_counts[class_counts >= TISSUE_MIN_PER_CLASS]
            eligible_core = len(strong) >= 2 and n_studies >= TISSUE_MIN_STUDIES_PER_CLASS
            if not eligible_core:
                reason = (
                    f"need ≥2 classes with ≥{TISSUE_MIN_PER_CLASS} samples and "
                    f"≥{TISSUE_MIN_STUDIES_PER_CLASS} studies"
                )
            eligible_aux = n_samples >= 100
        elif phenotype_id == "sex":
            eligible_aux = n_samples >= SEX_MIN_SAMPLES and n_studies >= 2
            eligible_core = False
            reason = "sex is auxiliary biological / QC, not a core burden target"
        elif phenotype_id == "ancestry":
            eligible_ext = n_samples >= 100
            reason = "ancestry is fairness / domain eval only"
        elif phenotype_id in {"blood", "brain"}:
            eligible_ext = n_samples >= 100
            eligible_aux = n_studies >= 2 and n_samples >= 200
            reason = "fine-grained blood/brain outside single-study default core"
        else:
            reason = "no eligibility rule for phenotype"

        rows.append(
            {
                "phenotype_id": phenotype_id,
                "phenotype_family": source_family,
                "task_type": task_type,
                "n_samples": n_samples,
                "n_cases": n_cases if phenotype_id in {"disease", "cancer"} else None,
                "n_controls": n_controls if phenotype_id in {"disease", "cancer"} else None,
                "n_unknown": n_unknown,
                "prevalence": prevalence,
                "n_studies": n_studies,
                "n_platforms": n_platforms,
                "n_tissues": len(tissues),
                "eligible_core_task": bool(eligible_core),
                "eligible_auxiliary_task": bool(eligible_aux),
                "eligible_external_evaluation": bool(eligible_ext),
                "exclusion_reason": reason,
            }
        )
    return pd.DataFrame(rows)


def head_training_allowed(
    eligibility: pd.DataFrame,
    *,
    phenotype_id: str,
    phenotype_family: str,
) -> bool:
    """True when catalog allows training this head (core or auxiliary)."""
    row = eligibility.loc[
        (eligibility["phenotype_id"] == phenotype_id)
        & (eligibility["phenotype_family"] == phenotype_family)
    ]
    if row.empty:
        return False
    record = row.iloc[0]
    return bool(record["eligible_core_task"]) or bool(record["eligible_auxiliary_task"])


def head_ranking_eligible(
    eligibility: pd.DataFrame,
    *,
    phenotype_id: str,
    phenotype_family: str,
) -> bool:
    """True only for core-eligible heads (architecture ranking tables)."""
    row = eligibility.loc[
        (eligibility["phenotype_id"] == phenotype_id)
        & (eligibility["phenotype_family"] == phenotype_family)
    ]
    if row.empty:
        return False
    return bool(row.iloc[0]["eligible_core_task"])


def validate_multitask_head_eligibility(
    *,
    data_root: Path,
    disease_enabled: bool,
    cancer_enabled: bool,
    disease_family: str = "disease",
    cancer_family: str = "cancer",
    skip_check: bool = False,
) -> dict[str, Any]:
    """Validate disease/cancer heads against trait_eligibility; return metrics metadata."""
    if skip_check or (not disease_enabled and not cancer_enabled):
        return {}
    path = release_paths(data_root).catalog_tables / "trait_eligibility.parquet"
    if not path.is_file():
        raise FileNotFoundError(
            f"trait_eligibility missing at {path}; run mbs catalog refresh-release"
        )
    elig = pd.read_parquet(path)
    meta: dict[str, Any] = {}
    if disease_enabled:
        if not head_training_allowed(elig, phenotype_id="disease", phenotype_family=disease_family):
            raise ValueError(
                f"disease head not eligible for family {disease_family!r}; "
                "see reports/inspection/deepmat_data_v1/trait_eligibility.md"
            )
        meta["disease"] = {
            "phenotype_family": disease_family,
            "training_allowed": True,
            "ranking_eligible": head_ranking_eligible(
                elig, phenotype_id="disease", phenotype_family=disease_family
            ),
        }
    if cancer_enabled:
        if not head_training_allowed(elig, phenotype_id="cancer", phenotype_family=cancer_family):
            raise ValueError(
                f"cancer head not eligible for family {cancer_family!r}; "
                "see reports/inspection/deepmat_data_v1/trait_eligibility.md"
            )
        meta["cancer"] = {
            "phenotype_family": cancer_family,
            "training_allowed": True,
            "ranking_eligible": head_ranking_eligible(
                elig, phenotype_id="cancer", phenotype_family=cancer_family
            ),
        }
    return meta


def validate_release_manifest(manifest: dict[str, Any]) -> None:
    required = [
        "artifact_version",
        "release_id",
        "genome_build",
        "created_at",
        "git_commit",
        "phenotype_families",
        "dedup_policy",
        "preprocessing",
        "probe_universe",
        "graph_id",
        "static_feature_set_id",
        "source_checksums",
        "ewas_db",
    ]
    missing = [key for key in required if key not in manifest]
    if missing:
        raise ValueError(f"release manifest missing keys: {missing}")
    if manifest["artifact_version"] != ARTIFACT_VERSION:
        raise ValueError(f"artifact_version must be {ARTIFACT_VERSION}")
    if manifest["genome_build"] != "GRCh38":
        raise ValueError("genome_build must be GRCh38")
    if not _COMMIT_RE.fullmatch(str(manifest["git_commit"])):
        raise ValueError("git_commit must be a 40-char lowercase hex SHA")
    if not isinstance(manifest["phenotype_families"], list) or not manifest["phenotype_families"]:
        raise ValueError("phenotype_families must be a non-empty list")
    if not isinstance(manifest["source_checksums"], list):
        raise TypeError("source_checksums must be a list")
    for entry in manifest["source_checksums"]:
        if "path" not in entry or "sha256" not in entry:
            raise ValueError("source_checksums entry requires path and sha256")
        if not _SHA256_RE.fullmatch(str(entry["sha256"])):
            raise ValueError(f"invalid sha256 for {entry.get('path')}")
    ewas = manifest["ewas_db"]
    for key in (
        "n_local_studies",
        "n_local_gsm",
        "advertised_n",
        "mirror_complete",
        "refresh_is_idempotent",
    ):
        if key not in ewas:
            raise ValueError(f"ewas_db missing {key}")


def _write_parquet(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(path, index=False)


def _load_table(connection: duckdb.DuckDBPyConnection, table: str, path: Path) -> None:
    if not path.is_file():
        return
    frame = pd.read_parquet(path)
    if frame.empty:
        return
    connection.execute(f"DELETE FROM {table}")  # noqa: S608
    connection.register("_mbs_load", frame)
    try:
        connection.execute(
            f"INSERT INTO {table} BY NAME SELECT * FROM _mbs_load"  # noqa: S608
        )
    finally:
        connection.unregister("_mbs_load")


def _empty_with_columns(columns: list[str]) -> pd.DataFrame:
    return pd.DataFrame({col: pd.Series(dtype="object") for col in columns})


def _populate_duckdb(
    *,
    database: Path,
    sql_dir: Path,
    tables_dir: Path,
) -> None:
    if database.exists():
        database.unlink()
    build_catalog(
        database=database,
        sql_dir=sql_dir,
        parquet_root=tables_dir,
        read_only=False,
    )
    connection = duckdb.connect(str(database))
    try:
        order = [
            "source_release",
            "platform",
            "study",
            "sample",
            "assay_file",
            "phenotype",
            "sample_phenotype",
            "sample_source_membership",
            "matrix_artifact",
            "matrix_sample",
            "fold_assignment",
            "artifact",
            "experiment",
            "trait_eligibility",
            "study_atlas_enrichment",
        ]
        for table in order:
            _load_table(connection, table, tables_dir / f"{table}.parquet")
    finally:
        connection.close()


def refresh_release(
    *,
    paths: DataPaths,
    release_id: str = RELEASE_ID,
    fetch_remote_index: bool = False,
    report_dir: Path | None = None,
) -> RefreshResult:
    """Build or refresh ``deepmat-data-v1`` from Hub Parquet + EWAS_db listing."""
    paths.ensure_directories()
    (paths.data_root / "canonical" / "releases").mkdir(parents=True, exist_ok=True)
    rp = release_paths(paths.data_root, release_id)
    for directory in (
        rp.root,
        rp.catalog_tables,
        rp.phenotypes_dir,
        rp.matrices_dir,
        rp.ontologies_dir,
        rp.splits_dir,
    ):
        directory.mkdir(parents=True, exist_ok=True)

    phenotypes, membership, hub_samples = build_long_form_phenotypes(paths.data_root)
    ewas_files = scan_ewas_db_tree(paths.data_root / "raw" / "ewas_datahub" / "EWAS_db")
    study_inv = ewas_db_study_inventory(ewas_files)
    matrix_artifacts, matrix_samples = scan_matrix_artifacts(paths.data_root)
    fold_rows = load_frozen_split(paths.artifact_root)

    advertised_n = EWAS_DB_ADVERTISED_N
    remote_names: list[str] | None = None
    if fetch_remote_index:
        remote_names = fetch_remote_ewas_db_study_names()
        if remote_names:
            advertised_n = len(remote_names)
            _write_parquet(
                rp.catalog_tables / "ewas_db_remote_studies.parquet",
                pd.DataFrame({"study_id": remote_names}),
            )

    # Studies / samples union
    platforms = _seed_platforms()
    known_platforms = set(platforms["platform_id"]) | set(PLATFORM_ALIASES.values())

    source_releases = pd.DataFrame(
        [
            {
                "source_release_id": HUB_SOURCE_RELEASE_ID,
                "source_name": "EWAS DataHub Baseline packs",
                "source_system": "ewas_datahub_baseline",
                "source_version": "v1",
                "retrieved_at": utc_now_iso(),
                "source_uri": "https://download.cncb.ac.cn/ewas/datahub/download/",
                "license_note": "EWAS Open Platform public download",
                "manifest_sha256": None,
            },
            {
                "source_release_id": EWAS_DB_SOURCE_RELEASE_ID,
                "source_name": "EWAS DataHub All Data",
                "source_system": "ewas_datahub_db",
                "source_version": "mirror-in-progress",
                "retrieved_at": utc_now_iso(),
                "source_uri": EWAS_DB_INDEX_URL,
                "license_note": "EWAS Open Platform public download",
                "manifest_sha256": None,
            },
        ]
    )

    study_rows: dict[str, dict[str, Any]] = {}
    if not hub_samples.empty:
        for rec in hub_samples.to_dict(orient="records"):
            study_id = str(rec["study_id"])
            platform_id = rec.get("platform_id")
            if platform_id is not None and platform_id not in known_platforms:
                platform_id = None
            study_rows[study_id] = {
                "study_id": study_id,
                "source_release_id": HUB_SOURCE_RELEASE_ID,
                "gse_id": study_id,
                "cohort_id": None,
                "platform_id": platform_id,
                "processing_level": "GMQN",
                "genome_build": "GRCh38",
                "retrieved_at": utc_now_iso(),
                "metadata_json": json.dumps({"lanes": ["ewas_datahub_baseline"]}),
            }
    for study_id in study_inv["study_id"].tolist() if not study_inv.empty else []:
        if study_id in study_rows:
            meta = json.loads(study_rows[study_id]["metadata_json"])
            lanes = list(dict.fromkeys([*meta.get("lanes", []), "ewas_datahub_db"]))
            study_rows[study_id]["metadata_json"] = json.dumps({"lanes": lanes})
            continue
        study_rows[study_id] = {
            "study_id": study_id,
            "source_release_id": EWAS_DB_SOURCE_RELEASE_ID,
            "gse_id": study_id,
            "cohort_id": None,
            "platform_id": None,
            "processing_level": "raw_beta_txt",
            "genome_build": "GRCh38",
            "retrieved_at": utc_now_iso(),
            "metadata_json": json.dumps({"lanes": ["ewas_datahub_db"]}),
        }
    studies = pd.DataFrame(list(study_rows.values())) if study_rows else pd.DataFrame()

    atlas_enrichment = build_study_atlas_enrichment(
        catalog_study_ids=[str(s) for s in studies["study_id"].tolist()] if not studies.empty else [],
        atlas_root=paths.data_root / "raw" / "ewas_atlas",
        gse_map_path=paths.project_root / "configs/data/atlas_gse_es_map.tsv",
    )
    if not studies.empty and not atlas_enrichment.empty:
        studies = merge_atlas_enrichment_into_studies(studies, atlas_enrichment)

    sample_frames: list[pd.DataFrame] = []
    if not hub_samples.empty:
        hub_sample_frame = pd.DataFrame(
            [
                {
                    "sample_id": str(rec["sample_id"]),
                    "study_id": str(rec["study_id"]),
                    "source_sample_id": str(rec.get("source_sample_id") or rec["sample_id"]),
                    "donor_id": rec.get("donor_id"),
                    "replicate_group": rec.get("replicate_group"),
                    "age": rec.get("age"),
                    "sex": rec.get("sex"),
                    "tissue_raw": rec.get("tissue_raw"),
                    "tissue_ontology_id": None,
                    "case_control": rec.get("case_control"),
                    "metadata_json": None,
                }
                for rec in hub_samples.to_dict(orient="records")
            ]
        )
        sample_frames.append(hub_sample_frame)
    if not ewas_files.empty:
        hub_ids = set(sample_frames[0]["sample_id"].astype(str)) if sample_frames else set()
        ewas_only = ewas_files.loc[
            ~ewas_files["sample_id"].astype(str).isin(list(hub_ids)),
            ["sample_id", "study_id"],
        ].drop_duplicates(subset=["sample_id"])
        if not ewas_only.empty:
            sid_series = ewas_only["sample_id"].astype(str)
            study_series = ewas_only["study_id"].astype(str)
            sample_frames.append(
                pd.DataFrame(
                    {
                        "sample_id": sid_series.to_numpy(),
                        "study_id": study_series.to_numpy(),
                        "source_sample_id": sid_series.to_numpy(),
                        "donor_id": None,
                        "replicate_group": None,
                        "age": None,
                        "sex": None,
                        "tissue_raw": None,
                        "tissue_ontology_id": None,
                        "case_control": None,
                        "metadata_json": json.dumps({"source": "ewas_db"}),
                    }
                )
            )
    known_sample_ids: set[str] = set()
    for frame in sample_frames:
        known_sample_ids.update(frame["sample_id"].astype(str).tolist())
    # Ensure matrix samples exist
    if not matrix_samples.empty:
        matrix_only_rows: list[dict[str, Any]] = []
        for rec in matrix_samples.to_dict(orient="records"):
            sid = str(rec["sample_id"])
            if sid in known_sample_ids:
                continue
            known_sample_ids.add(sid)
            sentinel = "__matrix_only__"
            if sentinel not in study_rows:
                study_rows[sentinel] = {
                    "study_id": sentinel,
                    "source_release_id": HUB_SOURCE_RELEASE_ID,
                    "gse_id": None,
                    "cohort_id": None,
                    "platform_id": None,
                    "processing_level": "matrix_pointer",
                    "genome_build": "GRCh38",
                    "retrieved_at": utc_now_iso(),
                    "metadata_json": json.dumps({"lanes": ["matrix_index"]}),
                }
                studies = pd.DataFrame(list(study_rows.values()))
            matrix_only_rows.append(
                {
                    "sample_id": sid,
                    "study_id": sentinel,
                    "source_sample_id": str(rec.get("source_sample_id") or sid),
                    "donor_id": None,
                    "replicate_group": None,
                    "age": None,
                    "sex": None,
                    "tissue_raw": None,
                    "tissue_ontology_id": None,
                    "case_control": None,
                    "metadata_json": json.dumps({"source": "matrix_index"}),
                }
            )
        if matrix_only_rows:
            sample_frames.append(pd.DataFrame(matrix_only_rows))
    samples = pd.concat(sample_frames, ignore_index=True) if sample_frames else pd.DataFrame()
    if not samples.empty:
        samples = samples.drop_duplicates(subset=["sample_id"], keep="first")

    if not ewas_files.empty:
        assay_files = pd.DataFrame(
            {
                "assay_file_id": (
                    "ewasdb:"
                    + ewas_files["study_id"].astype(str)
                    + ":"
                    + ewas_files["sample_id"].astype(str)
                ),
                "study_id": ewas_files["study_id"].astype(str),
                "path": ewas_files["path"].astype(str),
                "format": "probe_beta_txt",
                "sha256": ewas_files["sha256"].astype(str),
                "byte_size": ewas_files["byte_size"].astype("int64"),
                "n_rows": pd.Series([None] * len(ewas_files), dtype="object"),
                "n_columns": 2,
                "matrix_orientation": pd.Series([None] * len(ewas_files), dtype="object"),
                "schema_hash": pd.Series([None] * len(ewas_files), dtype="object"),
                "processing_level": "raw_beta_txt",
            }
        )
    else:
        assay_files = pd.DataFrame()

    phenotype_dim = _seed_phenotype_table()
    # Drop phenotype observations whose sample vanished (should not happen)
    if not phenotypes.empty and not samples.empty:
        phenotypes = _as_dataframe(phenotypes[phenotypes["sample_id"].isin(samples["sample_id"])])
    if not membership.empty and not samples.empty:
        membership = _as_dataframe(membership[membership["sample_id"].isin(samples["sample_id"])])
    if not matrix_samples.empty and not samples.empty:
        matrix_samples = _as_dataframe(
            matrix_samples[matrix_samples["sample_id"].isin(samples["sample_id"])]
        )
    if not fold_rows.empty and not samples.empty:
        fold_rows = _as_dataframe(fold_rows[fold_rows["sample_id"].isin(samples["sample_id"])])

    geo_merge_stats: dict[str, Any] = {"enabled": False}
    geo_path = paths.data_root / "canonical" / "phenotypes" / GEO_PARQUET_NAME
    if geo_backfill_enabled() and geo_path.is_file():
        geo_frame = load_geo_frame(paths.data_root)
        samples, phenotypes, studies, geo_merge_stats = merge_geo_sample_metadata(
            samples=samples,
            phenotypes=phenotypes,
            studies=studies,
            geo_frame=geo_frame,
        )
        geo_merge_stats["enabled"] = True
        geo_merge_stats["parquet_path"] = str(geo_path)

    if not samples.empty and study_rows:
        samples_for_elig = samples.copy()
        platform_lookup = {
            str(sid): (meta or {}).get("platform_id") for sid, meta in study_rows.items()
        }
        samples_for_elig["platform_id"] = samples_for_elig["study_id"].map(
            lambda sid: platform_lookup.get(str(sid))  # type: ignore[misc]
        )
        eligibility = compute_trait_eligibility(phenotypes, samples_for_elig)
    else:
        eligibility = compute_trait_eligibility(phenotypes, samples)

    # Artifacts / experiment pointers
    artifacts = pd.DataFrame(
        [
            {
                "artifact_id": f"release:{release_id}",
                "artifact_type": "data_release",
                "path": str(rp.root),
                "sha256": None,
                "manifest_json": json.dumps({"release_id": release_id}),
                "created_at": utc_now_iso(),
            }
        ]
    )
    experiments = pd.DataFrame(
        [
            {
                "experiment_id": f"catalog-refresh:{release_id}",
                "git_commit": git_commit(paths.project_root),
                "resolved_config_path": str(
                    paths.project_root / "docs" / "plans" / "milestone-7a-harmonized-release.md"
                ),
                "data_release_id": release_id,
                "graph_artifact_id": None,
                "feature_artifact_ids": [],
                "split_id": (
                    str(pd.Series(fold_rows["split_id"]).iloc[0]) if not fold_rows.empty else None
                ),
                "status": "complete",
                "created_at": utc_now_iso(),
            }
        ]
    )

    def _ensure(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
        if frame.empty:
            return _empty_with_columns(columns)
        return frame

    studies = _ensure(
        studies,
        [
            "study_id",
            "source_release_id",
            "gse_id",
            "cohort_id",
            "platform_id",
            "processing_level",
            "genome_build",
            "retrieved_at",
            "metadata_json",
        ],
    )
    samples = _ensure(
        samples,
        [
            "sample_id",
            "study_id",
            "source_sample_id",
            "donor_id",
            "replicate_group",
            "age",
            "sex",
            "tissue_raw",
            "tissue_ontology_id",
            "case_control",
            "metadata_json",
        ],
    )
    assay_files = _ensure(
        assay_files,
        [
            "assay_file_id",
            "study_id",
            "path",
            "format",
            "sha256",
            "byte_size",
            "n_rows",
            "n_columns",
            "matrix_orientation",
            "schema_hash",
            "processing_level",
        ],
    )
    phenotypes = _ensure(
        phenotypes,
        [
            "sample_id",
            "phenotype_id",
            "numeric_value",
            "categorical_value",
            "label_status",
            "is_observed",
            "source_family",
            "source_record_id",
            "ontology_id",
        ],
    )
    membership = _ensure(
        membership,
        [
            "sample_id",
            "source_release_id",
            "phenotype_family",
            "source_file",
            "source_row",
            "matrix_id",
            "row_index",
        ],
    )
    matrix_artifacts = _ensure(
        matrix_artifacts,
        [
            "matrix_id",
            "path",
            "platform_id",
            "processing_level",
            "genome_build",
            "n_samples",
            "n_loci",
            "manifest_sha256",
            "phenotype_family",
            "notes",
        ],
    )
    matrix_samples = _ensure(
        matrix_samples,
        ["matrix_id", "sample_id", "row_index", "source_sample_id"],
    )
    fold_rows = _ensure(
        fold_rows,
        ["split_id", "sample_id", "outer_fold", "inner_fold", "group_id", "role"],
    )
    eligibility = _ensure(
        eligibility,
        [
            "phenotype_id",
            "phenotype_family",
            "task_type",
            "n_samples",
            "n_cases",
            "n_controls",
            "n_unknown",
            "prevalence",
            "n_studies",
            "n_platforms",
            "n_tissues",
            "eligible_core_task",
            "eligible_auxiliary_task",
            "eligible_external_evaluation",
            "exclusion_reason",
        ],
    )
    atlas_enrichment = _ensure(
        atlas_enrichment,
        [
            "study_id",
            "join_method",
            "atlas_study_ids",
            "pmid",
            "n_atlas_cohorts",
            "total_sample_size",
            "tissues",
            "cohort_descriptions",
            "platforms",
            "ancestries",
            "atlas_traits",
        ],
    )

    # Durable parquet
    _write_parquet(rp.catalog_tables / "source_release.parquet", source_releases)
    _write_parquet(rp.catalog_tables / "platform.parquet", platforms)
    _write_parquet(rp.catalog_tables / "study.parquet", studies)
    _write_parquet(rp.catalog_tables / "sample.parquet", samples)
    _write_parquet(rp.catalog_tables / "assay_file.parquet", assay_files)
    _write_parquet(rp.catalog_tables / "phenotype.parquet", phenotype_dim)
    _write_parquet(rp.catalog_tables / "sample_phenotype.parquet", phenotypes)
    _write_parquet(rp.catalog_tables / "sample_source_membership.parquet", membership)
    _write_parquet(rp.catalog_tables / "matrix_artifact.parquet", matrix_artifacts)
    _write_parquet(rp.catalog_tables / "matrix_sample.parquet", matrix_samples)
    _write_parquet(rp.catalog_tables / "fold_assignment.parquet", fold_rows)
    _write_parquet(rp.catalog_tables / "artifact.parquet", artifacts)
    _write_parquet(rp.catalog_tables / "experiment.parquet", experiments)
    _write_parquet(rp.catalog_tables / "trait_eligibility.parquet", eligibility)
    _write_parquet(rp.catalog_tables / "ewas_db_study_inventory.parquet", study_inv)
    _write_parquet(rp.catalog_tables / "study_atlas_enrichment.parquet", atlas_enrichment)
    _write_parquet(rp.phenotypes_dir / "sample_phenotype.parquet", phenotypes)
    _write_parquet(rp.phenotypes_dir / "sample_source_membership.parquet", membership)
    _write_parquet(rp.matrices_dir / "index.parquet", matrix_artifacts)

    # Ontologies / split copies (small)
    for name in (
        "tissue_ontology.yaml",
        "tissue_ontology_age_tissue_sex_full_v1.yaml",
        "sex_ontology_v1.yaml",
    ):
        src = paths.data_root / "canonical" / "phenotypes" / name
        if src.is_file():
            shutil.copy2(src, rp.ontologies_dir / name)
    split_src = paths.artifact_root / "runs" / FROZEN_SPLIT_RUN_ID / "split.json"
    if split_src.is_file():
        dest = rp.splits_dir / "stage0-flat-deeprvat-age-tissue-sex-full-v1.split.json"
        shutil.copy2(split_src, dest)

    sql_dir = paths.project_root / "sql"
    _populate_duckdb(database=rp.catalog_db, sql_dir=sql_dir, tables_dir=rp.catalog_tables)

    n_local_studies = int(pd.Series(study_inv["study_id"]).nunique()) if not study_inv.empty else 0
    n_local_gsm = len(ewas_files)
    total_bytes = int(pd.Series(ewas_files["byte_size"]).sum()) if not ewas_files.empty else 0
    mirror_complete = n_local_studies >= advertised_n > 0

    source_checksums: list[dict[str, Any]] = []
    for family in HUB_FAMILIES:
        path = paths.data_root / "canonical" / "phenotypes" / f"{family}_sample_info.parquet"
        if path.is_file():
            source_checksums.append(
                {
                    "path": str(path),
                    "sha256": sha256_file(path),
                    "role": f"hub_sample_info:{family}",
                    "byte_size": int(path.stat().st_size),
                }
            )
    hub_dl = paths.data_root / "raw" / "ewas_datahub" / "download"
    if hub_dl.is_dir():
        for zip_path in sorted(hub_dl.glob("*_methylation_v1.zip")):
            size = int(zip_path.stat().st_size)
            source_checksums.append(
                {
                    "path": str(zip_path.resolve()),
                    "sha256": path_size_digest(zip_path.name, size),
                    "role": "hub_profile_or_sample_zip",
                    "byte_size": size,
                    "sha256_note": "content-address of name+size (pack too large to hash)",
                }
            )
    graph_manifest = (
        paths.data_root / "canonical" / "graphs" / DEFAULT_GRAPH_ID / "graph_manifest.json"
    )
    if graph_manifest.is_file():
        source_checksums.append(
            {
                "path": str(graph_manifest),
                "sha256": sha256_file(graph_manifest),
                "role": "annotation_graph_manifest",
                "byte_size": int(graph_manifest.stat().st_size),
            }
        )
    static_artifact = (
        paths.data_root / "canonical" / "static_features" / DEFAULT_FEATURE_SET_ID / "artifact.json"
    )
    if static_artifact.is_file():
        source_checksums.append(
            {
                "path": str(static_artifact),
                "sha256": sha256_file(static_artifact),
                "role": "static_feature_artifact",
                "byte_size": int(static_artifact.stat().st_size),
            }
        )
    inv_path = rp.catalog_tables / "ewas_db_study_inventory.parquet"
    if inv_path.is_file():
        source_checksums.append(
            {
                "path": str(inv_path),
                "sha256": sha256_file(inv_path),
                "role": "ewas_db_study_inventory",
                "byte_size": int(inv_path.stat().st_size),
            }
        )
    if geo_path.is_file():
        source_checksums.append(
            {
                "path": str(geo_path),
                "sha256": sha256_file(geo_path),
                "role": "geo_sample_metadata_backfill",
                "byte_size": int(geo_path.stat().st_size),
            }
        )

    present_families = sorted(
        set(membership["phenotype_family"].tolist()) if not membership.empty else []
    )
    manifest: dict[str, Any] = {
        "artifact_version": ARTIFACT_VERSION,
        "release_id": release_id,
        "genome_build": "GRCh38",
        "created_at": utc_now_iso(),
        "git_commit": git_commit(paths.project_root),
        "phenotype_families": present_families or list(HUB_FAMILIES),
        "dedup_policy": (
            "unique GSM sample_id; long-form sample_phenotype keeps multi-pack rows "
            "keyed by (sample_id, phenotype_id, source_family)"
        ),
        "preprocessing": {
            "hub_baseline_packs": "GMQN (DataHub published)",
            "ewas_db": "raw probe_id\\tbeta text (no clipping)",
        },
        "probe_universe": DEFAULT_GRAPH_ID,
        "graph_id": DEFAULT_GRAPH_ID,
        "static_feature_set_id": DEFAULT_FEATURE_SET_ID,
        "source_checksums": source_checksums,
        "ewas_db": {
            "n_local_studies": n_local_studies,
            "n_local_gsm": n_local_gsm,
            "advertised_n": advertised_n,
            "mirror_complete": mirror_complete,
            "refresh_is_idempotent": True,
            "total_bytes": total_bytes,
            "remote_index_fetched": remote_names is not None,
        },
        "geo_backfill": geo_merge_stats,
        "catalog_path": str(rp.catalog_db),
        "notes": (
            "Re-run mbs catalog refresh-release after EWAS_db download adds study dirs. "
            "Does not overwrite deepmat-data-age-tissue-sex-v1 or v0.1 model freezes."
        ),
    }
    validate_release_manifest(manifest)
    write_json(rp.manifest_path, manifest)

    resolved_report = report_dir
    census_before = load_census_snapshot(resolved_report) if resolved_report is not None else None
    if resolved_report is not None:
        write_phenotype_census_report(
            database=rp.catalog_db,
            report_dir=resolved_report,
            release_manifest=manifest,
        )
        write_trait_eligibility_report(
            database=rp.catalog_db,
            report_dir=resolved_report,
        )
        write_study_atlas_enrichment_report(
            enrichment=atlas_enrichment,
            report_dir=resolved_report,
        )
        if geo_merge_stats.get("enabled"):
            write_geo_backfill_pilot_report(
                stats=geo_merge_stats,
                database=rp.catalog_db,
                report_dir=resolved_report,
                census_before=census_before,
                census_after=load_census_snapshot(resolved_report),
            )

    return RefreshResult(
        release_id=release_id,
        release_root=str(rp.root),
        catalog_path=str(rp.catalog_db),
        n_samples=len(samples),
        n_studies=len(studies),
        n_phenotype_rows=len(phenotypes),
        ewas_db_n_local_studies=n_local_studies,
        ewas_db_n_local_gsm=n_local_gsm,
        ewas_db_mirror_complete=mirror_complete,
        report_dir=str(resolved_report) if resolved_report is not None else None,
    )


def validate_release(*, data_root: Path, release_id: str = RELEASE_ID) -> dict[str, Any]:
    rp = release_paths(data_root, release_id)
    if not rp.manifest_path.is_file():
        raise FileNotFoundError(f"release manifest missing: {rp.manifest_path}")
    if not rp.catalog_db.is_file():
        raise FileNotFoundError(f"release catalog missing: {rp.catalog_db}")
    manifest = json.loads(rp.manifest_path.read_text(encoding="utf-8"))
    validate_release_manifest(manifest)
    connection = duckdb.connect(str(rp.catalog_db), read_only=True)
    try:
        sample_row = connection.execute("SELECT count(*) FROM sample").fetchone()
        pheno_row = connection.execute("SELECT count(*) FROM sample_phenotype").fetchone()
        elig_row = connection.execute("SELECT count(*) FROM trait_eligibility").fetchone()
        if sample_row is None or pheno_row is None or elig_row is None:
            raise RuntimeError("catalog count queries returned no rows")
        n_samples = int(sample_row[0])
        n_pheno = int(pheno_row[0])
        n_elig = int(elig_row[0])
    finally:
        connection.close()
    return {
        "release_id": release_id,
        "manifest_ok": True,
        "catalog_path": str(rp.catalog_db),
        "n_samples": n_samples,
        "n_phenotype_rows": n_pheno,
        "n_trait_eligibility_rows": n_elig,
        "ewas_db": manifest.get("ewas_db"),
    }


def _query_df(database: Path, sql: str) -> pd.DataFrame:
    connection = duckdb.connect(str(database), read_only=True)
    try:
        return connection.execute(sql).fetchdf()
    finally:
        connection.close()


def write_phenotype_census_report(
    *,
    database: Path,
    report_dir: Path,
    release_manifest: dict[str, Any] | None = None,
) -> Path:
    report_dir = report_dir.resolve()
    report_dir.mkdir(parents=True, exist_ok=True)
    unique_gsm = int(
        _query_df(database, "SELECT count(DISTINCT sample_id) AS n FROM sample")["n"].iloc[0]
    )
    pack_rows = _query_df(
        database,
        """
        SELECT phenotype_family, count(*) AS n_rows,
               count(DISTINCT sample_id) AS n_unique_gsm
        FROM sample_source_membership
        GROUP BY phenotype_family
        ORDER BY phenotype_family
        """,
    )
    pack_row_sum = int(pd.Series(pack_rows["n_rows"]).sum()) if not pack_rows.empty else 0
    overlap = _query_df(
        database,
        """
        SELECT n_families, count(*) AS n_samples
        FROM v_sample_pack_overlap
        GROUP BY n_families
        ORDER BY n_families
        """,
    )
    conflicts = _query_df(database, "SELECT * FROM v_sample_label_conflicts LIMIT 50")
    prevalence = _query_df(database, "SELECT * FROM v_phenotype_prevalence ORDER BY 1, 2")
    ewas = _query_df(database, "SELECT * FROM v_ewas_db_ingest_status")
    age_by_study = _query_df(
        database,
        "SELECT * FROM v_age_distribution_by_study ORDER BY study_id LIMIT 500",
    )
    bmi_by_study = _query_df(
        database,
        "SELECT * FROM v_bmi_distribution_by_study ORDER BY study_id LIMIT 500",
    )
    donor_counts = _query_df(
        database,
        """
        SELECT
            count(*) FILTER (WHERE donor_id IS NOT NULL) AS n_with_donor,
            count(*) FILTER (WHERE replicate_group IS NOT NULL) AS n_with_replicate,
            count(DISTINCT donor_id) FILTER (WHERE donor_id IS NOT NULL) AS n_distinct_donors,
            count(DISTINCT replicate_group) FILTER (
                WHERE replicate_group IS NOT NULL
            ) AS n_distinct_replicate_groups
        FROM sample
        """,
    )
    payload = {
        "generated_at": utc_now_iso(),
        "unique_gsm": unique_gsm,
        "pack_row_sum": pack_row_sum,
        "unique_vs_row_sum_note": (
            "unique GSM must be compared to pack-row sum; overlap makes them unequal"
        ),
        "pack_membership": pack_rows.to_dict(orient="records"),
        "overlap_by_n_families": overlap.to_dict(orient="records"),
        "label_conflicts_head": conflicts.to_dict(orient="records"),
        "phenotype_prevalence": prevalence.to_dict(orient="records"),
        "ewas_db_ingest": ewas.to_dict(orient="records"),
        "release_ewas_db": (release_manifest or {}).get("ewas_db"),
        "age_distribution_by_study": age_by_study.to_dict(orient="records"),
        "bmi_distribution_by_study": bmi_by_study.to_dict(orient="records"),
        "donor_replicate": (
            donor_counts.to_dict(orient="records")[0] if not donor_counts.empty else {}
        ),
    }
    write_json(report_dir / "census.json", payload)
    lines = [
        "# Phenotype census (deepmat-data-v1)",
        "",
        f"- Generated: `{payload['generated_at']}`",
        f"- Unique GSM (`sample`): **{unique_gsm}**",
        f"- Pack membership row sum: **{pack_row_sum}**",
        "- Pack row sum counts Hub membership only; unique GSM also includes "
        "EWAS_db-only samples. Pack row sum can exceed unique Hub GSMs when "
        "samples appear in multiple packs.",
        "",
        "## Pack membership",
        "",
        "| Family | Rows | Unique GSM |",
        "| --- | ---: | ---: |",
    ]
    lines.extend(
        f"| `{row['phenotype_family']}` | {row['n_rows']} | {row['n_unique_gsm']} |"
        for row in cast(list[dict[str, Any]], payload["pack_membership"])
    )
    lines.extend(["", "## Overlap by number of packs", ""])
    if overlap.empty:
        lines.append("_No membership rows._")
    else:
        lines.append("| n_families | n_samples |")
        lines.append("| ---: | ---: |")
        lines.extend(
            f"| {row['n_families']} | {row['n_samples']} |"
            for row in cast(list[dict[str, Any]], payload["overlap_by_n_families"])
        )
    lines.extend(
        [
            "",
            "## Label conflicts (head)",
            "",
            f"Rows: {len(conflicts)} (capped at 50 in report).",
            "",
            "## EWAS_db ingest",
            "",
        ]
    )
    ewas_meta = cast(dict[str, Any], payload.get("release_ewas_db") or {})
    lines.append(
        f"- Local studies: `{ewas_meta.get('n_local_studies')}` / "
        f"advertised `{ewas_meta.get('advertised_n')}` "
        f"(mirror_complete={ewas_meta.get('mirror_complete')})"
    )
    lines.append(f"- Local GSM files: `{ewas_meta.get('n_local_gsm')}`")
    lines.append("")
    donor_meta = cast(dict[str, Any], payload.get("donor_replicate") or {})
    lines.extend(
        [
            "## Donor / replicate (when present)",
            "",
            f"- Samples with donor_id: **{donor_meta.get('n_with_donor', 0)}**",
            f"- Distinct donors: **{donor_meta.get('n_distinct_donors', 0)}**",
            f"- Samples with replicate_group: **{donor_meta.get('n_with_replicate', 0)}**",
            f"- Distinct replicate groups: **{donor_meta.get('n_distinct_replicate_groups', 0)}**",
            "",
            "Hub sample-info currently lacks donor columns for most packs; "
            "counts stay 0 until source fields exist (never invented from GSM).",
            "",
            "## Within-study age ranges",
            "",
        ]
    )
    age_rows = cast(list[dict[str, Any]], payload.get("age_distribution_by_study") or [])
    if not age_rows:
        lines.append("_No age-by-study rows._")
    else:
        lines.append("| study_id | n | age_min | age_max | age_mean |")
        lines.append("| --- | ---: | ---: | ---: | ---: |")
        lines.extend(
            f"| `{row['study_id']}` | {row['n_samples']} | "
            f"{row['age_min']} | {row['age_max']} | {row['age_mean']} |"
            for row in age_rows[:50]
        )
        if len(age_rows) > 50:
            lines.append(f"| … | ({len(age_rows) - 50} more studies) | | | |")
    lines.extend(["", "## Within-study BMI ranges", ""])
    bmi_rows = cast(list[dict[str, Any]], payload.get("bmi_distribution_by_study") or [])
    if not bmi_rows:
        lines.append("_No BMI-by-study rows._")
    else:
        lines.append("| study_id | n | bmi_min | bmi_max | bmi_mean |")
        lines.append("| --- | ---: | ---: | ---: | ---: |")
        lines.extend(
            f"| `{row['study_id']}` | {row['n_samples']} | "
            f"{row['bmi_min']} | {row['bmi_max']} | {row['bmi_mean']} |"
            for row in bmi_rows[:50]
        )
        if len(bmi_rows) > 50:
            lines.append(f"| … | ({len(bmi_rows) - 50} more studies) | | | |")
    lines.append("")
    lines.append("Re-run `mbs catalog refresh-release` after more `EWAS_db` study dirs download.")
    lines.append("")
    (report_dir / "census.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report_dir


def write_trait_eligibility_report(*, database: Path, report_dir: Path) -> Path:
    report_dir = report_dir.resolve()
    report_dir.mkdir(parents=True, exist_ok=True)
    frame = _query_df(database, "SELECT * FROM v_trait_training_eligibility")
    rows = cast(list[dict[str, Any]], frame.to_dict(orient="records"))
    payload = {
        "generated_at": utc_now_iso(),
        "rows": rows,
    }
    write_json(report_dir / "trait_eligibility.json", payload)
    lines = [
        "# Trait eligibility (deepmat-data-v1)",
        "",
        f"- Generated: `{payload['generated_at']}`",
        "",
        "| Family | Phenotype | Task | n | studies | core | aux | ext | Reason |",
        "| --- | --- | --- | ---: | ---: | --- | --- | --- | --- |",
    ]
    lines.extend(
        (
            "| `{phenotype_family}` | `{phenotype_id}` | {task_type} | {n_samples} | "
            "{n_studies} | {eligible_core_task} | {eligible_auxiliary_task} | "
            "{eligible_external_evaluation} | {exclusion_reason} |"
        ).format(
            phenotype_family=row["phenotype_family"],
            phenotype_id=row["phenotype_id"],
            task_type=row["task_type"],
            n_samples=row["n_samples"],
            n_studies=row["n_studies"],
            eligible_core_task=row["eligible_core_task"],
            eligible_auxiliary_task=row["eligible_auxiliary_task"],
            eligible_external_evaluation=row["eligible_external_evaluation"],
            exclusion_reason=(
                "" if _is_blank(row.get("exclusion_reason")) else str(row["exclusion_reason"])
            ),
        )
        for row in rows
    )
    lines.append("")
    lines.append(
        "Disease/cancer: unknown labels are never treated as controls (ADR / DATA_CONTRACT)."
    )
    lines.append("")
    (report_dir / "trait_eligibility.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report_dir
