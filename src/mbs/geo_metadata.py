"""Sample-level GEO SOFT metadata backfill for EWAS_db-only GSM (not Atlas labels)."""

from __future__ import annotations

import gzip
import json
import os
import re
import time
import urllib.request
from pathlib import Path
from typing import Any, cast

import duckdb
import pandas as pd

from mbs.annotation.manifest import sha256_file, utc_now_iso, write_json

GEO_SOURCE_FAMILY = "geo_metadata_backfill"
GEO_PARQUET_NAME = "geo_sample_metadata.parquet"
GEO_MANIFEST_NAME = "geo_sample_metadata.manifest.json"
NCBI_AGENT = "methyl-burden-score/1.0 (research pipeline)"
NCBI_DELAY_S = 0.34
CATALOG_PLATFORMS = frozenset({"HM450", "EPIC", "EPICv2"})

# NCBI GEO methylation BeadChip accessions → catalog platform_id.
# Do not map expression / unknown GPLs (they stay in sample.metadata_json.geo.platform_id).
GPL_TO_PLATFORM: dict[str, str] = {
    "GPL13534": "HM450",  # Illumina HumanMethylation450
    "GPL16304": "HM450",  # 450K (UBC enhanced annotation; alt to GPL13534)
    "GPL21145": "EPIC",  # Infinium MethylationEPIC
    "GPL23976": "EPICv2",  # HumanMethylation850 / EPIC 850k
}

_AGE_SENTINELS = frozenset({10002, 999, 9999})
_CONTROL_TOKENS = frozenset(
    {"control", "controls", "healthy", "normal", "unaffected", "non-disease", "non disease"}
)
_CASE_TOKENS = frozenset({"case", "patient", "diseased", "affected", "tumor", "tumour"})
_CANCER_KEY_RE = re.compile(r"cancer|tumor|tumour|malignan|carcinoma|neoplasm", re.I)
_BATCH_KEY_RE = re.compile(r"^batch|batch id|batch_id|lot|plate", re.I)
_TREATMENT_KEY_RE = re.compile(r"^treatment|therapy|drug|medication", re.I)
_SOFT_TABLE_BEGIN = re.compile(r"^!(?:Sample_|Series_)?(?:platform|sample)_table_begin", re.I)
_SOFT_TABLE_END = re.compile(r"^!(?:Sample_|Series_)?(?:platform|sample)_table_end", re.I)


def family_soft_url(gse_id: str) -> str:
    """NCBI FTP URL for a series family SOFT file."""
    gse = gse_id.strip().upper()
    if not gse.startswith("GSE"):
        raise ValueError(f"not a GSE accession: {gse_id!r}")
    suffix = gse.removeprefix("GSE")
    if not suffix.isdigit():
        raise ValueError(f"not a GSE accession: {gse_id!r}")
    stem = suffix[:-3] if len(suffix) > 3 else "0"
    return (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE{stem}nnn/{gse}/soft/{gse}_family.soft.gz"
    )


def cache_soft_path(cache_root: Path, gse_id: str) -> Path:
    gse = gse_id.strip().upper()
    return cache_root / "geo" / gse / f"{gse}_family.soft.gz"


def geo_parquet_path(data_root: Path) -> Path:
    return data_root / "canonical" / "phenotypes" / GEO_PARQUET_NAME


def geo_manifest_path(data_root: Path) -> Path:
    return data_root / "canonical" / "phenotypes" / GEO_MANIFEST_NAME


def geo_backfill_enabled() -> bool:
    return os.environ.get("MBS_SKIP_GEO_BACKFILL", "0") != "1"


def catalog_platform_from_gpl(gpl: str | None) -> str | None:
    """Map a GEO GPL accession to catalog platform_id, or None if unknown."""
    if gpl is None or (isinstance(gpl, float) and pd.isna(gpl)):
        return None
    text = str(gpl).strip().upper()
    if not text:
        return None
    mapped = GPL_TO_PLATFORM.get(text)
    if mapped in CATALOG_PLATFORMS:
        return mapped
    return None


def _is_blank(value: object) -> bool:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return True
    return str(value).strip() == ""


def _parse_soft_blocks(text: str) -> tuple[dict[str, str], list[dict[str, list[str]]]]:
    """Parse family SOFT; skip embedded platform/sample tables."""
    series: dict[str, str] = {}
    samples: list[dict[str, list[str]]] = []
    current: dict[str, list[str]] | None = None
    in_table = False
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if _SOFT_TABLE_BEGIN.match(line):
            in_table = True
            continue
        if _SOFT_TABLE_END.match(line):
            in_table = False
            continue
        if in_table:
            continue
        if line.startswith("^SERIES"):
            current = None
            continue
        if line.startswith("^SAMPLE"):
            if current:
                samples.append(current)
            current = {}
            continue
        if not line.startswith("!"):
            continue
        body = line[1:].strip()
        if "=" not in body:
            continue
        key, value = body.split("=", 1)
        key = key.strip()
        value = value.strip()
        if current is None:
            series.setdefault(key, value)
        else:
            current.setdefault(key, []).append(value)
    if current:
        samples.append(current)
    return series, samples


def _normalize_key(key: str) -> str:
    text = key.strip().lower()
    for prefix in ("sample_", "series_"):
        if text.startswith(prefix):
            text = text[len(prefix) :]
    return text


def _parse_characteristics(sample: dict[str, list[str]]) -> dict[str, str]:
    chars: dict[str, str] = {}
    for key, values in sample.items():
        norm = _normalize_key(key)
        if norm in {"characteristics_ch1", "characteristics_ch2"} or norm.startswith(
            "characteristics_"
        ):
            for value in values:
                if ":" in value:
                    sub_key, sub_val = value.split(":", 1)
                    chars[sub_key.strip().lower()] = sub_val.strip()
                else:
                    chars[norm] = value.strip()
    return chars


def _parse_age(raw: str) -> tuple[float | None, str, bool]:
    text = raw.strip()
    if not text:
        return None, "unknown", False
    match = re.search(r"(\d+(?:\.\d+)?)", text)
    if not match:
        return None, "unknown", False
    age = float(match.group(1))
    if age < 0 or age > 120 or int(age) in _AGE_SENTINELS:
        return None, "unknown", False
    return age, "observed", True


def _parse_sex(raw: str) -> tuple[str | None, str, bool]:
    text = raw.strip()
    if not text:
        return None, "unknown", False
    lowered = text.lower()
    if lowered in {"m", "male", "man"}:
        return "Male", "observed", True
    if lowered in {"f", "female", "woman"}:
        return "Female", "observed", True
    return None, "unknown", False


def _disease_case_control(raw: str, *, key: str) -> tuple[str | None, str | None, bool]:
    """Return (categorical_value, label_status, write_row) for disease/cancer."""
    text = raw.strip()
    if not text:
        return None, None, False
    lowered = text.lower()
    tokens = set(re.split(r"[\s,;/]+", lowered))
    is_control = bool(tokens & _CONTROL_TOKENS) or lowered.startswith("control")
    is_case = bool(tokens & _CASE_TOKENS) or lowered in {"case", "patient", "affected"}
    if not is_control and not is_case:
        return None, None, False
    if is_control:
        return text, "control", True
    return text, "case", True


def characteristics_to_phenotypes(chars: dict[str, str]) -> list[dict[str, Any]]:
    """Map GEO characteristics to catalog phenotype observations."""
    rows: list[dict[str, Any]] = []
    for key, value in chars.items():
        if _BATCH_KEY_RE.search(key) or _TREATMENT_KEY_RE.search(key):
            continue
        norm = key.strip().lower()
        if norm in {"age", "age (years)", "age years"}:
            numeric, status, observed = _parse_age(value)
            if observed and numeric is not None:
                rows.append(
                    {
                        "phenotype_id": "age",
                        "numeric_value": numeric,
                        "categorical_value": None,
                        "label_status": status,
                        "is_observed": True,
                    }
                )
            continue
        if norm in {"sex", "gender"}:
            cat, status, observed = _parse_sex(value)
            rows.append(
                {
                    "phenotype_id": "sex",
                    "numeric_value": None,
                    "categorical_value": cat,
                    "label_status": status,
                    "is_observed": observed,
                }
            )
            continue
        if norm in {"tissue", "cell type", "cell_type", "organism part", "organ", "source tissue"}:
            if not _is_blank(value):
                rows.append(
                    {
                        "phenotype_id": "tissue",
                        "numeric_value": None,
                        "categorical_value": value.strip(),
                        "label_status": "observed",
                        "is_observed": True,
                    }
                )
            continue
        if norm in {
            "disease",
            "diagnosis",
            "condition",
            "disease status",
            "disease state",
            "health status",
        } or "disease" in norm or "diagnosis" in norm:
            cat, status, write = _disease_case_control(value, key=norm)
            if write and cat is not None and status is not None:
                phenotype_id = "cancer" if _CANCER_KEY_RE.search(norm) else "disease"
                rows.append(
                    {
                        "phenotype_id": phenotype_id,
                        "numeric_value": None,
                        "categorical_value": cat,
                        "label_status": status,
                        "is_observed": True,
                    }
                )
    return rows


def _sample_field(sample: dict[str, list[str]], *names: str) -> str | None:
    for name in names:
        for key, values in sample.items():
            if _normalize_key(key) == _normalize_key(name) and values:
                text = values[0].strip()
                if text:
                    return text
    return None


def parse_family_soft(text: str) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Parse family SOFT into series metadata and per-GSM sample records."""
    series_raw, sample_blocks = _parse_soft_blocks(text)
    series_id = series_raw.get("Series_geo_accession") or series_raw.get("geo_accession")
    pubmed_raw = series_raw.get("Series_pubmed_id") or series_raw.get("pubmed_id") or ""
    pubmed_ids = sorted(
        dict.fromkeys(p.strip() for p in re.split(r"[;,]\s*", pubmed_raw) if p.strip())
    )
    series = {
        "study_id": (series_id or "").strip().upper() or None,
        "pubmed_ids": pubmed_ids,
        "title": series_raw.get("Series_title") or series_raw.get("title"),
    }
    samples: list[dict[str, Any]] = []
    for block in sample_blocks:
        sample_id = _sample_field(block, "Sample_geo_accession", "geo_accession")
        if not sample_id:
            continue
        chars = _parse_characteristics(block)
        phenotypes = characteristics_to_phenotypes(chars)
        gpl = _sample_field(block, "Sample_platform_id", "platform_id")
        study_id_val = series["study_id"] or _sample_field(block, "Sample_series_id") or ""
        row: dict[str, Any] = {
            "sample_id": sample_id.strip().upper(),
            "study_id": study_id_val.upper(),
            "source_name": _sample_field(block, "Sample_source_name_ch1", "source_name_ch1"),
            "platform_id": gpl,
            "catalog_platform_id": catalog_platform_from_gpl(gpl),
            "pubmed_ids": pubmed_ids,
            "characteristics_raw": chars,
            "phenotypes": phenotypes,
        }
        for pheno in phenotypes:
            pid = pheno["phenotype_id"]
            if pid == "age" and pheno.get("is_observed"):
                row["age"] = pheno["numeric_value"]
            elif pid == "sex":
                row["sex"] = pheno.get("categorical_value")
                row["sex_label_status"] = pheno["label_status"]
            elif pid == "tissue" and pheno.get("is_observed"):
                row["tissue"] = pheno.get("categorical_value")
            elif pid in {"disease", "cancer"}:
                row[pid] = pheno.get("categorical_value")
                row[f"{pid}_label_status"] = pheno["label_status"]
        samples.append(row)
    return series, samples


def download_family_soft(
    gse_id: str,
    *,
    cache_root: Path,
    delay_s: float = NCBI_DELAY_S,
    force: bool = False,
) -> tuple[Path, str]:
    """Download and cache family SOFT; return path and sha256 of cached gzip."""
    cache_path = cache_soft_path(cache_root, gse_id)
    if cache_path.is_file() and not force:
        return cache_path, sha256_file(cache_path)
    url = family_soft_url(gse_id)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": NCBI_AGENT})  # noqa: S310
    with urllib.request.urlopen(req, timeout=120) as resp:  # noqa: S310
        data = resp.read()
    cache_path.write_bytes(data)
    if delay_s > 0:
        time.sleep(delay_s)
    return cache_path, sha256_file(cache_path)


def read_cached_soft(path: Path) -> str:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            return handle.read()
    return path.read_text(encoding="utf-8", errors="replace")


def build_geo_frame_from_soft(text: str, *, fetched_at: str, soft_sha256: str) -> pd.DataFrame:
    """Build one-row-per-GSM frame ready for Parquet export."""
    _, samples = parse_family_soft(text)
    rows: list[dict[str, Any]] = []
    for sample in samples:
        chars = sample.get("characteristics_raw") or {}
        rows.append(
            {
                "sample_id": sample["sample_id"],
                "study_id": sample["study_id"],
                "source_name": sample.get("source_name"),
                "platform_id": sample.get("platform_id"),
                "catalog_platform_id": sample.get("catalog_platform_id"),
                "pubmed_ids": json.dumps(sample.get("pubmed_ids") or []),
                "characteristics_raw": json.dumps(chars, sort_keys=True),
                "age": sample.get("age"),
                "sex": sample.get("sex"),
                "sex_label_status": sample.get("sex_label_status"),
                "tissue": sample.get("tissue"),
                "disease": sample.get("disease"),
                "disease_label_status": sample.get("disease_label_status"),
                "cancer": sample.get("cancer"),
                "cancer_label_status": sample.get("cancer_label_status"),
                "fetched_at": fetched_at,
                "soft_sha256": soft_sha256,
            }
        )
    return pd.DataFrame(rows)


def _hub_sample_ids(samples: pd.DataFrame) -> set[str]:
    """GSM present from Hub packs (metadata_json is null, not ewas_db-only)."""
    if samples.empty:
        return set()
    hub: set[str] = set()
    for rec in samples.to_dict(orient="records"):
        meta_raw = rec.get("metadata_json")
        if meta_raw is None or (isinstance(meta_raw, float) and pd.isna(meta_raw)):
            hub.add(str(rec["sample_id"]))
            continue
        if isinstance(meta_raw, str):
            try:
                meta = json.loads(meta_raw)
            except json.JSONDecodeError:
                continue
            if meta.get("source") != "ewas_db":
                hub.add(str(rec["sample_id"]))
    return hub


def _parse_metadata_json(raw: object) -> dict[str, Any]:
    if raw is None or (isinstance(raw, float) and pd.isna(raw)):
        return {}
    if isinstance(raw, dict):
        return dict(raw)
    if isinstance(raw, str) and raw.strip():
        return cast(dict[str, Any], json.loads(raw))
    return {}


def _geo_rows_to_phenotypes(row: dict[str, Any]) -> list[dict[str, Any]]:
    pheno_rows: list[dict[str, Any]] = []
    if row.get("age") is not None and not _is_blank(row.get("age")):
        pheno_rows.append(
            {
                "sample_id": row["sample_id"],
                "phenotype_id": "age",
                "numeric_value": float(row["age"]),
                "categorical_value": None,
                "label_status": "observed",
                "is_observed": True,
                "source_family": GEO_SOURCE_FAMILY,
                "source_record_id": f"geo:{row['sample_id']}:age",
                "ontology_id": None,
            }
        )
    if not _is_blank(row.get("sex")):
        status = str(row.get("sex_label_status") or "observed")
        pheno_rows.append(
            {
                "sample_id": row["sample_id"],
                "phenotype_id": "sex",
                "numeric_value": None,
                "categorical_value": str(row["sex"]),
                "label_status": status,
                "is_observed": status != "unknown",
                "source_family": GEO_SOURCE_FAMILY,
                "source_record_id": f"geo:{row['sample_id']}:sex",
                "ontology_id": None,
            }
        )
    if not _is_blank(row.get("tissue")):
        pheno_rows.append(
            {
                "sample_id": row["sample_id"],
                "phenotype_id": "tissue",
                "numeric_value": None,
                "categorical_value": str(row["tissue"]),
                "label_status": "observed",
                "is_observed": True,
                "source_family": GEO_SOURCE_FAMILY,
                "source_record_id": f"geo:{row['sample_id']}:tissue",
                "ontology_id": None,
            }
        )
    for pid in ("disease", "cancer"):
        if not _is_blank(row.get(pid)):
            status = str(row.get(f"{pid}_label_status") or "case")
            pheno_rows.append(
                {
                    "sample_id": row["sample_id"],
                    "phenotype_id": pid,
                    "numeric_value": None,
                    "categorical_value": str(row[pid]),
                    "label_status": status,
                    "is_observed": True,
                    "source_family": GEO_SOURCE_FAMILY,
                    "source_record_id": f"geo:{row['sample_id']}:{pid}",
                    "ontology_id": None,
                }
            )
    return pheno_rows


def merge_geo_sample_metadata(
    *,
    samples: pd.DataFrame,
    phenotypes: pd.DataFrame,
    studies: pd.DataFrame,
    geo_frame: pd.DataFrame,
    skip: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Merge GEO backfill rows for EWAS_db-only GSM; Hub GSM omitted entirely."""
    stats: dict[str, Any] = {
        "enabled": not skip,
        "n_geo_rows_input": len(geo_frame),
        "n_samples_touched": 0,
        "n_samples_skipped_hub": 0,
        "n_samples_skipped_missing": 0,
        "n_phenotype_rows_added": 0,
        "phenotypes_by_id": {},
        "studies_touched": [],
        "per_study": {},
    }
    if skip or geo_frame.empty:
        return samples, phenotypes, studies, stats

    hub_ids = _hub_sample_ids(samples)
    catalog_ids = set(samples["sample_id"].astype(str)) if not samples.empty else set()
    geo_frame = geo_frame.copy()
    geo_frame["sample_id"] = geo_frame["sample_id"].astype(str).str.strip().str.upper()
    if "platform_id" in geo_frame.columns:
        geo_frame["catalog_platform_id"] = [
            catalog_platform_from_gpl(v) for v in geo_frame["platform_id"].tolist()
        ]
    else:
        geo_frame["catalog_platform_id"] = None

    new_pheno_rows: list[dict[str, Any]] = []
    sample_index = (
        samples.set_index("sample_id", drop=False) if not samples.empty else pd.DataFrame()
    )
    studies = studies.copy() if not studies.empty else pd.DataFrame()
    study_index = studies.set_index("study_id", drop=False) if not studies.empty else None

    touched: set[str] = set()
    per_study: dict[str, dict[str, int]] = {}

    for rec in geo_frame.to_dict(orient="records"):
        sid = str(rec["sample_id"])
        if sid in hub_ids:
            stats["n_samples_skipped_hub"] += 1
            continue
        if sid not in catalog_ids:
            stats["n_samples_skipped_missing"] += 1
            continue
        pheno_rows = _geo_rows_to_phenotypes(rec)
        if not pheno_rows and _is_blank(rec.get("source_name")) and _is_blank(
            rec.get("characteristics_raw")
        ):
            continue
        touched.add(sid)
        study_id = str(rec.get("study_id") or "")
        per_study.setdefault(study_id, {"samples": 0, "phenotypes": 0})
        per_study[study_id]["samples"] += 1
        per_study[study_id]["phenotypes"] += len(pheno_rows)
        new_pheno_rows.extend(pheno_rows)
        for pheno in pheno_rows:
            pid = str(pheno["phenotype_id"])
            stats["phenotypes_by_id"][pid] = stats["phenotypes_by_id"].get(pid, 0) + 1

        if sid in sample_index.index:
            idx = sample_index.index.get_loc(sid)
            if isinstance(idx, slice):
                idx = idx.start or 0
            row = samples.loc[samples["sample_id"] == sid].iloc[0]
            if _is_blank(row.get("age")) and rec.get("age") is not None:
                samples.loc[samples["sample_id"] == sid, "age"] = rec["age"]
            if _is_blank(row.get("sex")) and not _is_blank(rec.get("sex")):
                samples.loc[samples["sample_id"] == sid, "sex"] = rec["sex"]
            if _is_blank(row.get("tissue_raw")) and not _is_blank(rec.get("tissue")):
                samples.loc[samples["sample_id"] == sid, "tissue_raw"] = rec["tissue"]
            for pid in ("disease", "cancer"):
                status_col = f"{pid}_label_status"
                if _is_blank(row.get("case_control")) and not _is_blank(rec.get(pid)):
                    cc = rec.get(status_col)
                    if cc in {"case", "control"}:
                        samples.loc[samples["sample_id"] == sid, "case_control"] = cc

            meta = _parse_metadata_json(row.get("metadata_json"))
            geo_meta = dict(meta.get("geo") or {})
            if not _is_blank(rec.get("source_name")):
                geo_meta["source_name"] = rec["source_name"]
            if not _is_blank(rec.get("platform_id")):
                geo_meta["platform_id"] = rec["platform_id"]
            if not _is_blank(rec.get("characteristics_raw")):
                try:
                    geo_meta["characteristics_raw"] = json.loads(str(rec["characteristics_raw"]))
                except json.JSONDecodeError:
                    geo_meta["characteristics_raw"] = rec["characteristics_raw"]
            geo_meta["fetched_at"] = rec.get("fetched_at")
            geo_meta["soft_sha256"] = rec.get("soft_sha256")
            meta["geo"] = geo_meta
            if meta.get("source") == "ewas_db" or "source" not in meta:
                meta["source"] = "ewas_db"
            samples.loc[samples["sample_id"] == sid, "metadata_json"] = json.dumps(
                meta, sort_keys=True
            )

        if study_index is not None and study_id and study_id in study_index.index:
            study_row = studies.loc[studies["study_id"] == study_id].iloc[0]
            study_meta = _parse_metadata_json(study_row.get("metadata_json"))
            geo_study = dict(study_meta.get("geo") or {})
            pubmed_ids: list[str] = []
            if not _is_blank(rec.get("pubmed_ids")):
                try:
                    pubmed_ids = json.loads(str(rec["pubmed_ids"]))
                except json.JSONDecodeError:
                    pubmed_ids = []
            if pubmed_ids:
                existing = geo_study.get("pubmed_ids") or []
                geo_study["pubmed_ids"] = sorted(dict.fromkeys([*existing, *pubmed_ids]))
                study_meta["geo"] = geo_study
                studies.loc[studies["study_id"] == study_id, "metadata_json"] = json.dumps(
                    study_meta, sort_keys=True
                )

    if study_index is not None and not geo_frame.empty:
        for study_id, grp in geo_frame.groupby(geo_frame["study_id"].astype(str)):
            if not study_id or study_id not in study_index.index:
                continue
            study_row = studies.loc[studies["study_id"] == study_id].iloc[0]
            if not _is_blank(study_row.get("platform_id")):
                continue
            plats = {
                p
                for p in grp["catalog_platform_id"].dropna().astype(str).tolist()
                if p in CATALOG_PLATFORMS
            }
            if len(plats) == 1:
                studies.loc[studies["study_id"] == study_id, "platform_id"] = next(iter(plats))

    stats["n_samples_touched"] = len(touched)
    stats["n_phenotype_rows_added"] = len(new_pheno_rows)
    stats["per_study"] = per_study
    stats["studies_touched"] = sorted(per_study.keys())

    if new_pheno_rows:
        new_phenotypes = pd.DataFrame(new_pheno_rows)
        phenotypes = (
            pd.concat([phenotypes, new_phenotypes], ignore_index=True)
            if not phenotypes.empty
            else new_phenotypes
        )
        phenotypes = phenotypes.drop_duplicates(
            subset=["sample_id", "phenotype_id", "source_family"],
            keep="first",
        )

    return samples, phenotypes, studies, stats


def load_geo_frame(data_root: Path) -> pd.DataFrame:
    path = geo_parquet_path(data_root)
    if not path.is_file():
        return pd.DataFrame()
    return pd.read_parquet(path)


def write_geo_parquet(data_root: Path, frame: pd.DataFrame, *, study_ids: list[str]) -> Path:
    out = geo_parquet_path(data_root)
    out.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(out, index=False)
    manifest = {
        "artifact_version": "geo-sample-metadata-v1",
        "created_at": utc_now_iso(),
        "n_samples": len(frame),
        "study_ids": sorted(set(study_ids)),
        "parquet_path": str(out),
        "sha256": sha256_file(out),
    }
    write_json(geo_manifest_path(data_root), manifest)
    return out


def compact_census_snapshot(census: dict[str, Any] | None) -> dict[str, Any] | None:
    """Keep census delta fields only (do not embed the full census blob)."""
    if not census:
        return None
    if "geo_prevalence" in census and "phenotype_prevalence" not in census:
        return census
    prevalence = census.get("phenotype_prevalence") or []
    geo_prev = [
        row
        for row in prevalence
        if isinstance(row, dict) and row.get("source_family") == GEO_SOURCE_FAMILY
    ]
    return {
        "generated_at": census.get("generated_at"),
        "unique_gsm": census.get("unique_gsm"),
        "pack_row_sum": census.get("pack_row_sum"),
        "geo_prevalence": geo_prev,
    }


def load_census_snapshot(report_dir: Path | None) -> dict[str, Any] | None:
    if report_dir is None:
        return None
    path = report_dir / "census.json"
    if not path.is_file():
        return None
    return compact_census_snapshot(json.loads(path.read_text(encoding="utf-8")))


def write_geo_backfill_pilot_report(
    *,
    stats: dict[str, Any],
    database: Path,
    report_dir: Path,
    census_before: dict[str, Any] | None = None,
    census_after: dict[str, Any] | None = None,
) -> Path:
    """Write pilot audit summary under geo_backfill_pilot/."""
    pilot_dir = report_dir / "geo_backfill_pilot"
    pilot_dir.mkdir(parents=True, exist_ok=True)
    connection = duckdb.connect(str(database), read_only=True)
    try:
        prevalence = connection.execute(
            """
            SELECT phenotype_id, source_family, n_rows, n_observed, n_samples
            FROM v_phenotype_prevalence
            WHERE source_family = ?
            ORDER BY phenotype_id
            """,
            [GEO_SOURCE_FAMILY],
        ).fetchdf()
        geo_label = connection.execute(
            """
            SELECT phenotype_id, label_status, count(*) AS n
            FROM sample_phenotype
            WHERE source_family = ?
            GROUP BY 1, 2
            ORDER BY 1, 2
            """,
            [GEO_SOURCE_FAMILY],
        ).fetchdf()
        geo_samples_with_pheno = connection.execute(
            """
            SELECT count(DISTINCT sample_id) AS n
            FROM sample_phenotype
            WHERE source_family = ? AND is_observed
            """,
            [GEO_SOURCE_FAMILY],
        ).fetchone()
        n_geo_observed = int(geo_samples_with_pheno[0]) if geo_samples_with_pheno else 0
        hub_overlap = connection.execute(
            """
            SELECT count(*) AS n
            FROM sample_phenotype sp
            JOIN sample_source_membership m USING (sample_id)
            WHERE sp.source_family = ?
            """,
            [GEO_SOURCE_FAMILY],
        ).fetchone()
        n_hub_overlap = int(hub_overlap[0]) if hub_overlap else 0
        n_geo_on_sample = connection.execute(
            """
            SELECT count(*) AS n
            FROM sample
            WHERE metadata_json IS NOT NULL
              AND json_extract(metadata_json, '$.geo') IS NOT NULL
            """
        ).fetchone()
        n_samples_with_geo_json = int(n_geo_on_sample[0]) if n_geo_on_sample else 0
        atlas_on_sample = connection.execute(
            """
            SELECT count(*) AS n
            FROM sample
            WHERE metadata_json IS NOT NULL
              AND json_extract(metadata_json, '$.atlas_enrichment') IS NOT NULL
            """
        ).fetchone()
        n_atlas_on_sample = int(atlas_on_sample[0]) if atlas_on_sample else 0
        eligibility = connection.execute(
            """
            SELECT phenotype_id, n_samples, n_cases, n_controls, n_unknown,
                   eligible_core_task, eligible_auxiliary_task, exclusion_reason
            FROM trait_eligibility
            WHERE phenotype_family = ?
            ORDER BY phenotype_id
            """,
            [GEO_SOURCE_FAMILY],
        ).fetchdf()
        membership_families = connection.execute(
            """
            SELECT phenotype_family, count(*) AS n
            FROM sample_source_membership
            GROUP BY 1
            ORDER BY 1
            """
        ).fetchdf()
    finally:
        connection.close()

    census_before = compact_census_snapshot(census_before) or census_before
    census_after = compact_census_snapshot(census_after) or census_after
    unique_before = (census_before or {}).get("unique_gsm")
    unique_after = (census_after or {}).get("unique_gsm")
    geo_prev_before = (census_before or {}).get("geo_prevalence") or []
    n_geo_obs_before = int(sum(int(r.get("n_observed") or 0) for r in geo_prev_before))

    summary: dict[str, Any] = {
        "generated_at": utc_now_iso(),
        "geo_source_family": GEO_SOURCE_FAMILY,
        "invariants": {
            "hub_wins_geo_rows_on_hub_membership_gsm": n_hub_overlap,
            "atlas_enrichment_on_sample": n_atlas_on_sample,
            "sample_source_membership_includes_geo": bool(
                not membership_families.empty
                and GEO_SOURCE_FAMILY in set(membership_families["phenotype_family"].astype(str))
            ),
        },
        "merge_stats": stats,
        "n_samples_with_geo_json": n_samples_with_geo_json,
        "n_samples_with_geo_phenotype": n_geo_observed,
        "prevalence": prevalence.to_dict(orient="records") if not prevalence.empty else [],
        "label_status": geo_label.to_dict(orient="records") if not geo_label.empty else [],
        "eligibility": eligibility.to_dict(orient="records") if not eligibility.empty else [],
        "census_before": census_before,
        "census_after": census_after,
        "census_delta": {
            "unique_gsm_before": unique_before,
            "unique_gsm_after": unique_after,
            "geo_observed_rows_before": n_geo_obs_before,
            "geo_observed_rows_after": int(prevalence["n_observed"].sum())
            if not prevalence.empty
            else 0,
            "n_geo_observed_samples": n_geo_observed,
            "phenotypes_by_id": stats.get("phenotypes_by_id", {}),
        },
        "notes": [
            "GEO rows are omitted for Hub GSM (Hub wins).",
            "Atlas enrichment stays on study.metadata_json only.",
            "Disease/cancer rows need explicit case/control tokens.",
            "Diagnosis-only text stays in metadata_json.",
            "Training heads read Hub pack Parquet, not geo_metadata_backfill.",
        ],
    }
    json_path = pilot_dir / "summary.json"
    write_json(json_path, summary)

    lines = [
        "# GEO backfill pilot summary",
        "",
        f"- Generated: `{summary['generated_at']}`",
        f"- GEO parquet GSM in: **{stats.get('n_geo_rows_input', 0)}**",
        f"- Catalog samples touched (EWAS_db-only): **{stats.get('n_samples_touched', 0)}**",
        f"- Hub-skipped GSM: **{stats.get('n_samples_skipped_hub', 0)}**",
        f"- GEO GSM not in catalog: **{stats.get('n_samples_skipped_missing', 0)}**",
        f"- Phenotype rows added: **{stats.get('n_phenotype_rows_added', 0)}**",
        f"- Samples with ≥1 observed GEO phenotype: **{n_geo_observed}**",
        f"- Samples with `metadata_json.geo`: **{n_samples_with_geo_json}**",
        "",
        "## Invariants",
        "",
        f"- GEO phenotype rows on Hub pack-membership GSM: **{n_hub_overlap}** (must be 0)",
        f"- Atlas blobs on `sample.metadata_json`: **{n_atlas_on_sample}** (must be 0)",
        f"- `sample_source_membership` includes `geo_metadata_backfill`: "
        f"**{summary['invariants']['sample_source_membership_includes_geo']}** (must be false)",
        "",
        "## Phenotypes by id",
        "",
        "| phenotype_id | rows | observed | unique GSM |",
        "| --- | ---: | ---: | ---: |",
    ]
    if prevalence.empty:
        lines.append("| _(none)_ | 0 | 0 | 0 |")
    else:
        lines.extend(
            (
                f"| `{rec['phenotype_id']}` | {rec['n_rows']} | "
                f"{rec['n_observed']} | {rec['n_samples']} |"
            )
            for rec in prevalence.to_dict(orient="records")
        )
    lines.extend(["", "## Label status", ""])
    if geo_label.empty:
        lines.append("_(none)_")
    else:
        lines.append("| phenotype_id | label_status | n |")
        lines.append("| --- | --- | ---: |")
        lines.extend(
            f"| `{rec['phenotype_id']}` | `{rec['label_status']}` | {rec['n']} |"
            for rec in geo_label.to_dict(orient="records")
        )
    lines.extend(["", "## Eligibility (`source_family=geo_metadata_backfill`)", ""])
    if eligibility.empty:
        lines.append("_(none)_")
    else:
        lines.append("| phenotype_id | n | cases | controls | core | aux |")
        lines.append("| --- | ---: | ---: | ---: | --- | --- |")
        for rec in eligibility.to_dict(orient="records"):
            core = rec["eligible_core_task"]
            aux = rec["eligible_auxiliary_task"]
            lines.append(
                f"| `{rec['phenotype_id']}` | {rec['n_samples']} | "
                f"{rec['n_cases']} | {rec['n_controls']} | {core} | {aux} |"
            )
            reason = rec.get("exclusion_reason")
            if reason and not rec.get("eligible_core_task"):
                lines.append(f"  - `{rec['phenotype_id']}` not core: {reason}")
    lines.extend(["", "## Per study (merge)", ""])
    lines.append("| study_id | samples touched | phenotype rows |")
    lines.append("| --- | ---: | ---: |")
    lines.extend(
        f"| `{study_id}` | {counts.get('samples', 0)} | {counts.get('phenotypes', 0)} |"
        for study_id, counts in sorted((stats.get("per_study") or {}).items())
    )
    geo_after = summary["census_delta"]["geo_observed_rows_after"]
    lines.extend(
        [
            "",
            "## Census delta",
            "",
            f"- Unique GSM before: **{unique_before}**",
            f"- Unique GSM after: **{unique_after}**",
            f"- GEO observed phenotype rows before: **{n_geo_obs_before}**",
            f"- GEO observed phenotype rows after: **{geo_after}**",
            "",
            "GEO backfill does not add `sample` rows (EWAS_db scan does). Unique-GSM "
            "movement is EWAS_db mirror growth, not GEO. Phenotype-row movement is the GEO delta.",
            "",
            "## Operator notes",
            "",
        ]
    )
    lines.extend(f"- {note}" for note in summary["notes"])
    md_path = pilot_dir / "summary.md"
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return pilot_dir


def load_census_before(project_root: Path) -> dict[str, Any] | None:
    """Deprecated path helper; prefer load_census_snapshot(report_dir)."""
    return load_census_snapshot(project_root / "reports" / "inspection" / "deepmat_data_v1")
