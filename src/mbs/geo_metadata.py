"""Sample-level GEO SOFT metadata backfill for EWAS_db-only GSM (not Atlas labels)."""

from __future__ import annotations

import gzip
import json
import os
import re
import time
import urllib.request
from pathlib import Path
from typing import Any, Protocol, cast

import duckdb
import pandas as pd
import yaml

from mbs.annotation.manifest import sha256_file, utc_now_iso, write_json

GEO_SOURCE_FAMILY = "geo_metadata_backfill"
GEO_PARQUET_NAME = "geo_sample_metadata.parquet"
GEO_MANIFEST_NAME = "geo_sample_metadata.manifest.json"
NCBI_AGENT = "methyl-burden-score/1.0 (research pipeline)"
NCBI_DELAY_S = 0.34
CATALOG_PLATFORMS = frozenset({"HM450", "EPIC", "EPICv2"})
DEFAULT_TISSUE_ONTOLOGY_REL = Path("canonical/phenotypes/tissue_ontology_hub_nine_pack_v1.yaml")
DEFAULT_TISSUE_ALIASES_REL = Path("configs/data/geo_tissue_aliases.yaml")

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
_AGE_UNIT_RE = re.compile(
    r"(\d+(?:\.\d+)?)\s*(years?|yrs?|y|months?|mos?|mo|weeks?|wks?|w|days?|d)\b",
    re.I,
)
_AGE_BARE_RE = re.compile(r"^(\d+(?:\.\d+)?)\s*$")
_UNIT_TO_YEARS: dict[str, float] = {
    "y": 1.0,
    "yr": 1.0,
    "yrs": 1.0,
    "year": 1.0,
    "years": 1.0,
    "mo": 1.0 / 12.0,
    "mos": 1.0 / 12.0,
    "month": 1.0 / 12.0,
    "months": 1.0 / 12.0,
    "w": 1.0 / 52.0,
    "wk": 1.0 / 52.0,
    "wks": 1.0 / 52.0,
    "week": 1.0 / 52.0,
    "weeks": 1.0 / 52.0,
    "d": 1.0 / 365.25,
    "day": 1.0 / 365.25,
    "days": 1.0 / 365.25,
}
_CONFLICT_FIELDS = (
    "age",
    "sex",
    "tissue",
    "tissue_ontology_id",
    "disease",
    "cancer",
    "platform_id",
    "catalog_platform_id",
)


class TissueOntologyLike(Protocol):
    """Minimal tissue ontology surface (avoids importing training → release)."""

    labels: tuple[str, ...]
    label_to_id: dict[str, int]


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


def _parse_age(raw: str) -> tuple[float | None, str, bool, str, str | None]:
    """Parse age to years; return (years, status, observed, age_raw, unit_tag).

    Unit tag is ``years`` / ``months`` / ``weeks`` / ``days`` / ``mixed`` / None.
    Bare numbers are treated as years. Original string is always returned as age_raw.
    """
    text = raw.strip()
    if not text:
        return None, "unknown", False, text, None
    unit_hits = _AGE_UNIT_RE.findall(text)
    years = 0.0
    units_seen: list[str] = []
    if unit_hits:
        for value_s, unit_s in unit_hits:
            factor = _UNIT_TO_YEARS.get(unit_s.lower())
            if factor is None:
                return None, "unknown", False, text, None
            years += float(value_s) * factor
            units_seen.append(unit_s.lower())
        unit_tag = units_seen[0] if len(set(units_seen)) == 1 else "mixed"
        # Canonical unit family for reporting.
        if unit_tag != "mixed":
            unit_tag = next(
                (
                    name
                    for name, keys in (
                        ("years", {"y", "yr", "yrs", "year", "years"}),
                        ("months", {"mo", "mos", "month", "months"}),
                        ("weeks", {"w", "wk", "wks", "week", "weeks"}),
                        ("days", {"d", "day", "days"}),
                    )
                    if unit_tag in keys
                ),
                unit_tag,
            )
    else:
        bare = _AGE_BARE_RE.match(text)
        if bare is None:
            # Legacy-ish: first number only when no unit token — reject unit-like junk.
            if re.search(r"[A-Za-z]", text):
                return None, "unknown", False, text, None
            match = re.search(r"(\d+(?:\.\d+)?)", text)
            if not match:
                return None, "unknown", False, text, None
            years = float(match.group(1))
            unit_tag = "years"
        else:
            years = float(bare.group(1))
            unit_tag = "years"
    if years < 0 or years > 120 or int(years) in _AGE_SENTINELS:
        return None, "unknown", False, text, unit_tag
    return years, "observed", True, text, unit_tag


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


def load_geo_tissue_aliases(path: Path | None = None) -> dict[str, str]:
    """Load raw→canonical tissue aliases (keys lowercased)."""
    alias_path = path
    if alias_path is None:
        root = Path(os.environ.get("MBS_ROOT", ".")).resolve()
        alias_path = root / DEFAULT_TISSUE_ALIASES_REL
    if not alias_path.is_file():
        return {}
    data = yaml.safe_load(alias_path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise TypeError(f"geo tissue aliases must be a mapping: {alias_path}")
    raw_aliases = data.get("aliases") or {}
    if not isinstance(raw_aliases, dict):
        raise TypeError(f"geo tissue aliases.aliases must be a mapping: {alias_path}")
    return {
        str(k).strip().lower(): str(v).strip()
        for k, v in raw_aliases.items()
        if str(k).strip()
    }


def resolve_tissue_ontology_path(data_root: Path | None = None) -> Path | None:
    """Prefer Hub nine-pack ontology, then ATS, then pilot yaml."""
    roots: list[Path] = []
    if data_root is not None:
        roots.append(data_root)
    env_root = os.environ.get("MBS_DATA_ROOT")
    if env_root:
        roots.append(Path(env_root))
    candidates = (
        "tissue_ontology_hub_nine_pack_v1.yaml",
        "tissue_ontology_age_tissue_sex_full_v1.yaml",
        "tissue_ontology.yaml",
    )
    for root in roots:
        for name in candidates:
            path = root / "canonical" / "phenotypes" / name
            if path.is_file():
                return path
    return None


def map_geo_tissue(
    raw: str,
    *,
    ontology: TissueOntologyLike | None,
    aliases: dict[str, str] | None = None,
) -> tuple[str | None, str | None, str]:
    """Map a GEO tissue string.

    Returns ``(canonical_label, ontology_id, status)`` where status is
    ``mapped`` / ``unmapped`` / ``ambiguous`` / ``empty``.
    """
    text = raw.strip()
    if not text:
        return None, None, "empty"
    if ontology is None:
        return text, None, "unmapped"
    alias_map = aliases if aliases is not None else {}
    key = text.lower()
    candidates: list[str] = []
    if key in alias_map:
        candidates.append(alias_map[key])
    label_by_fold = {lab.casefold(): lab for lab in ontology.labels}
    if key.casefold() in label_by_fold:
        candidates.append(label_by_fold[key.casefold()])
    # Deduplicate while preserving order.
    uniq = list(dict.fromkeys(candidates))
    if len(uniq) > 1:
        return None, None, "ambiguous"
    if not uniq:
        return None, None, "unmapped"
    label = uniq[0]
    if label not in ontology.label_to_id:
        # Alias pointed at a label absent from this ontology.
        fold_hit = label_by_fold.get(label.casefold())
        if fold_hit is None:
            return None, None, "unmapped"
        label = fold_hit
    return label, str(ontology.label_to_id[label]), "mapped"

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
        if norm in {"age", "age (years)", "age years", "age (months)", "age months"}:
            numeric, status, observed, age_raw, age_unit = _parse_age(value)
            if observed and numeric is not None:
                rows.append(
                    {
                        "phenotype_id": "age",
                        "numeric_value": numeric,
                        "categorical_value": None,
                        "label_status": status,
                        "is_observed": True,
                        "age_raw": age_raw,
                        "age_unit": age_unit,
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
                        "tissue_raw": value.strip(),
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
                row["age_raw"] = pheno.get("age_raw")
                row["age_unit"] = pheno.get("age_unit")
            elif pid == "sex":
                row["sex"] = pheno.get("categorical_value")
                row["sex_label_status"] = pheno["label_status"]
            elif pid == "tissue" and pheno.get("is_observed"):
                row["tissue_raw"] = pheno.get("tissue_raw") or pheno.get("categorical_value")
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


def build_geo_frame_from_soft(
    text: str,
    *,
    fetched_at: str,
    soft_sha256: str,
    ontology: TissueOntologyLike | None = None,
    aliases: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Build one-row-per-GSM frame ready for Parquet export."""
    _, samples = parse_family_soft(text)
    alias_map = aliases if aliases is not None else load_geo_tissue_aliases()
    tissue_stats = {"mapped": 0, "unmapped": 0, "ambiguous": 0, "empty": 0}
    rows: list[dict[str, Any]] = []
    for sample in samples:
        chars = sample.get("characteristics_raw") or {}
        tissue_raw = sample.get("tissue_raw") or sample.get("tissue")
        tissue_label = None
        tissue_ont = None
        tissue_status = "empty"
        if not _is_blank(tissue_raw):
            tissue_label, tissue_ont, tissue_status = map_geo_tissue(
                str(tissue_raw),
                ontology=ontology,
                aliases=alias_map,
            )
            tissue_stats[tissue_status] = tissue_stats.get(tissue_status, 0) + 1
            if tissue_status != "mapped":
                tissue_label = None
                tissue_ont = None
        rows.append(
            {
                "sample_id": sample["sample_id"],
                "study_id": sample["study_id"],
                "study_ids": json.dumps([sample["study_id"]] if sample.get("study_id") else []),
                "source_name": sample.get("source_name"),
                "platform_id": sample.get("platform_id"),
                "catalog_platform_id": sample.get("catalog_platform_id"),
                "pubmed_ids": json.dumps(sample.get("pubmed_ids") or []),
                "characteristics_raw": json.dumps(chars, sort_keys=True),
                "age": sample.get("age"),
                "age_raw": sample.get("age_raw"),
                "age_unit": sample.get("age_unit"),
                "sex": sample.get("sex"),
                "sex_label_status": sample.get("sex_label_status"),
                "tissue_raw": tissue_raw,
                "tissue": tissue_label if tissue_status == "mapped" else tissue_raw,
                "tissue_ontology_id": tissue_ont,
                "tissue_map_status": tissue_status,
                "disease": sample.get("disease"),
                "disease_label_status": sample.get("disease_label_status"),
                "cancer": sample.get("cancer"),
                "cancer_label_status": sample.get("cancer_label_status"),
                "fetched_at": fetched_at,
                "soft_sha256": soft_sha256,
            }
        )
    frame = pd.DataFrame(rows)
    frame.attrs["tissue_map_stats"] = tissue_stats
    return frame


def _conflict_norm(field: str, value: object) -> str:
    if field == "age":
        return f"{float(value):.6g}"
    text = str(value).strip()
    if field in {"sex", "tissue", "disease", "cancer"}:
        return text.casefold()
    return text


def consolidate_geo_sample_rows(frame: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    """One row per GSM; persist multi-GSE membership; null conflicting phenotypes.

    Does **not** silently ``keep=\"first\"`` when metadata disagrees.
    """
    stats: dict[str, Any] = {
        "n_input_rows": len(frame),
        "n_unique_gsm": 0,
        "n_multi_study_gsm": 0,
        "n_conflict_samples": 0,
        "n_conflict_fields": 0,
        "conflicts": [],
    }
    if frame.empty:
        return frame.copy(), stats

    out_rows: list[dict[str, Any]] = []
    conflict_examples: list[dict[str, Any]] = []
    for sample_id, grp in frame.groupby(frame["sample_id"].astype(str), sort=False):
        study_ids: list[str] = []
        for raw in grp["study_id"].tolist():
            if _is_blank(raw):
                continue
            sid = str(raw).strip().upper()
            if sid and sid not in study_ids:
                study_ids.append(sid)
        # Also fold any pre-existing study_ids JSON columns.
        if "study_ids" in grp.columns:
            for raw in grp["study_ids"].tolist():
                if _is_blank(raw):
                    continue
                try:
                    parsed = json.loads(str(raw)) if isinstance(raw, str) else raw
                except json.JSONDecodeError:
                    parsed = []
                if isinstance(parsed, list):
                    for item in parsed:
                        text = str(item).strip().upper()
                        if text and text not in study_ids:
                            study_ids.append(text)
        study_ids = sorted(study_ids)
        if len(study_ids) > 1:
            stats["n_multi_study_gsm"] += 1
        base = grp.iloc[0].to_dict()
        base["sample_id"] = str(sample_id).strip().upper()
        base["study_id"] = study_ids[0] if study_ids else base.get("study_id")
        base["study_ids"] = json.dumps(study_ids)
        field_conflicts: list[dict[str, Any]] = []
        for field in _CONFLICT_FIELDS:
            if field not in grp.columns:
                continue
            non_blank = [
                v for v in grp[field].tolist() if not _is_blank(v)
            ]
            if not non_blank:
                base[field] = None
                continue
            norms = [_conflict_norm(field, v) for v in non_blank]
            uniq_norms = list(dict.fromkeys(norms))
            if len(uniq_norms) > 1:
                field_conflicts.append(
                    {
                        "field": field,
                        "values": [str(v) for v in dict.fromkeys(str(x) for x in non_blank)],
                    }
                )
                base[field] = None
                if field == "sex":
                    base["sex_label_status"] = None
                if field == "disease":
                    base["disease_label_status"] = None
                if field == "cancer":
                    base["cancer_label_status"] = None
                if field == "tissue":
                    base["tissue_map_status"] = "ambiguous"
                    base["tissue_ontology_id"] = None
            else:
                # Prefer first non-blank original value.
                base[field] = non_blank[0]
        # Prefer first non-blank age_raw / tissue_raw for provenance.
        for prov in ("age_raw", "age_unit", "tissue_raw", "sex_label_status"):
            if prov not in grp.columns:
                continue
            vals = [v for v in grp[prov].tolist() if not _is_blank(v)]
            if vals:
                base[prov] = vals[0]
        if field_conflicts:
            stats["n_conflict_samples"] += 1
            stats["n_conflict_fields"] += len(field_conflicts)
            if len(conflict_examples) < 50:
                conflict_examples.append(
                    {
                        "sample_id": base["sample_id"],
                        "study_ids": study_ids,
                        "fields": field_conflicts,
                    }
                )
        out_rows.append(base)

    stats["n_unique_gsm"] = len(out_rows)
    stats["conflicts"] = conflict_examples
    return pd.DataFrame(out_rows), stats


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
    tissue_status = str(row.get("tissue_map_status") or "")
    tissue_value = row.get("tissue")
    if tissue_status == "mapped" and not _is_blank(tissue_value):
        pheno_rows.append(
            {
                "sample_id": row["sample_id"],
                "phenotype_id": "tissue",
                "numeric_value": None,
                "categorical_value": str(tissue_value),
                "label_status": "observed",
                "is_observed": True,
                "source_family": GEO_SOURCE_FAMILY,
                "source_record_id": f"geo:{row['sample_id']}:tissue",
                "ontology_id": (
                    str(row["tissue_ontology_id"])
                    if not _is_blank(row.get("tissue_ontology_id"))
                    else None
                ),
            }
        )
    elif tissue_status in {"", "unmapped"} and not _is_blank(
        row.get("tissue_raw") or tissue_value
    ):
        # Persist raw for audit/census; ontology_id stays null.
        pheno_rows.append(
            {
                "sample_id": row["sample_id"],
                "phenotype_id": "tissue",
                "numeric_value": None,
                "categorical_value": str(row.get("tissue_raw") or tissue_value),
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
    ontology: TissueOntologyLike | None = None,
    aliases: dict[str, str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Merge GEO backfill rows for EWAS_db-only GSM; Hub GSM omitted entirely."""
    stats: dict[str, Any] = {
        "enabled": not skip,
        "n_geo_rows_input": len(geo_frame),
        "n_geo_phenotype_rows_before_merge": (
            int((phenotypes["source_family"] == GEO_SOURCE_FAMILY).sum())
            if not phenotypes.empty and "source_family" in phenotypes.columns
            else 0
        ),
        "n_samples_touched": 0,
        "n_samples_skipped_hub": 0,
        "n_samples_skipped_missing": 0,
        "n_phenotype_rows_added": 0,
        "phenotypes_by_id": {},
        "studies_touched": [],
        "per_study": {},
        "tissue_map": {"mapped": 0, "unmapped": 0, "ambiguous": 0, "empty": 0},
        "unmapped_tissue_examples": [],
        "n_multi_study_gsm": 0,
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

    alias_map = aliases if aliases is not None else load_geo_tissue_aliases()
    if ontology is not None and (
        "tissue_map_status" not in geo_frame.columns
        or geo_frame["tissue_map_status"].isna().all()
    ):
        mapped_tissue: list[object] = []
        mapped_ont: list[object] = []
        mapped_status: list[str] = []
        mapped_raw: list[object] = []
        for rec in geo_frame.to_dict(orient="records"):
            raw = rec.get("tissue_raw") or rec.get("tissue")
            mapped_raw.append(raw)
            if _is_blank(raw):
                mapped_tissue.append(None)
                mapped_ont.append(None)
                mapped_status.append("empty")
                continue
            label, ont_id, status = map_geo_tissue(
                str(raw), ontology=ontology, aliases=alias_map
            )
            mapped_status.append(status)
            if status == "mapped":
                mapped_tissue.append(label)
                mapped_ont.append(ont_id)
            else:
                mapped_tissue.append(None if status == "ambiguous" else raw)
                mapped_ont.append(None)
        geo_frame["tissue_raw"] = mapped_raw
        geo_frame["tissue"] = mapped_tissue
        geo_frame["tissue_ontology_id"] = mapped_ont
        geo_frame["tissue_map_status"] = mapped_status

    new_pheno_rows: list[dict[str, Any]] = []
    sample_index = (
        samples.set_index("sample_id", drop=False) if not samples.empty else pd.DataFrame()
    )
    studies = studies.copy() if not studies.empty else pd.DataFrame()
    study_index = studies.set_index("study_id", drop=False) if not studies.empty else None

    touched: set[str] = set()
    per_study: dict[str, dict[str, int]] = {}
    unmapped_examples: list[str] = []

    for rec in geo_frame.to_dict(orient="records"):
        sid = str(rec["sample_id"])
        if sid in hub_ids:
            stats["n_samples_skipped_hub"] += 1
            continue
        if sid not in catalog_ids:
            stats["n_samples_skipped_missing"] += 1
            continue
        t_status = str(rec.get("tissue_map_status") or "empty")
        if t_status in stats["tissue_map"]:
            stats["tissue_map"][t_status] += 1
        if t_status == "unmapped" and not _is_blank(rec.get("tissue_raw") or rec.get("tissue")):
            example = str(rec.get("tissue_raw") or rec.get("tissue"))
            if example not in unmapped_examples and len(unmapped_examples) < 40:
                unmapped_examples.append(example)
        study_ids_raw = rec.get("study_ids")
        study_id_list: list[str] = []
        if not _is_blank(study_ids_raw):
            try:
                parsed_ids = (
                    json.loads(str(study_ids_raw))
                    if isinstance(study_ids_raw, str)
                    else study_ids_raw
                )
            except json.JSONDecodeError:
                parsed_ids = []
            if isinstance(parsed_ids, list):
                study_id_list = [str(x).upper() for x in parsed_ids if not _is_blank(x)]
        if len(study_id_list) > 1:
            stats["n_multi_study_gsm"] += 1
        pheno_rows = _geo_rows_to_phenotypes(rec)
        if not pheno_rows and _is_blank(rec.get("source_name")) and _is_blank(
            rec.get("characteristics_raw")
        ):
            continue
        touched.add(sid)
        study_id = str(rec.get("study_id") or (study_id_list[0] if study_id_list else ""))
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
            tissue_for_sample = (
                rec.get("tissue")
                if str(rec.get("tissue_map_status") or "") == "mapped"
                else rec.get("tissue_raw") or rec.get("tissue")
            )
            if _is_blank(row.get("tissue_raw")) and not _is_blank(
                rec.get("tissue_raw") or tissue_for_sample
            ):
                samples.loc[samples["sample_id"] == sid, "tissue_raw"] = (
                    rec.get("tissue_raw") or tissue_for_sample
                )
            if _is_blank(row.get("tissue_ontology_id")) and not _is_blank(
                rec.get("tissue_ontology_id")
            ):
                samples.loc[samples["sample_id"] == sid, "tissue_ontology_id"] = rec[
                    "tissue_ontology_id"
                ]
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
            if not _is_blank(rec.get("age_raw")):
                geo_meta["age_raw"] = rec["age_raw"]
            if not _is_blank(rec.get("age_unit")):
                geo_meta["age_unit"] = rec["age_unit"]
            if not _is_blank(rec.get("tissue_raw")):
                geo_meta["tissue_raw"] = rec["tissue_raw"]
            if not _is_blank(rec.get("tissue_map_status")):
                geo_meta["tissue_map_status"] = rec["tissue_map_status"]
            if study_id_list:
                geo_meta["study_ids"] = study_id_list
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
    stats["unmapped_tissue_examples"] = unmapped_examples

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
    n_geo_obs_before_disk = int(sum(int(r.get("n_observed") or 0) for r in geo_prev_before))
    n_geo_before_merge = int(stats.get("n_geo_phenotype_rows_before_merge") or 0)
    dirty_disk_baseline = n_geo_obs_before_disk > 0 and n_geo_before_merge == 0

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
            "geo_phenotype_rows_before_merge": n_geo_before_merge,
            "geo_observed_rows_before_disk_census": n_geo_obs_before_disk,
            "geo_observed_rows_after": int(prevalence["n_observed"].sum())
            if not prevalence.empty
            else 0,
            "n_geo_observed_samples": n_geo_observed,
            "phenotypes_by_id": stats.get("phenotypes_by_id", {}),
            "dirty_disk_census_baseline": dirty_disk_baseline,
            "authoritative_delta_rows": int(stats.get("n_phenotype_rows_added") or 0),
        },
        "notes": [
            "GEO rows are omitted for Hub GSM (Hub wins).",
            "Atlas enrichment stays on study.metadata_json only.",
            "Disease/cancer rows need explicit case/control tokens.",
            "Diagnosis-only text stays in metadata_json.",
            "Training heads read Hub pack Parquet, not geo_metadata_backfill.",
            "Authoritative GEO Δ is merge_stats.n_phenotype_rows_added "
            "(in-memory phenotypes before merge are usually 0 on a full refresh).",
            "For a clean incremental test: MBS_SKIP_GEO_BACKFILL=1 refresh → assert "
            "zero geo_metadata_backfill rows → fetch/merge → compare exact Δ.",
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
        f"- Multi-study GSM (membership): **{stats.get('n_multi_study_gsm', 0)}**",
        "",
        "## Invariants",
        "",
        f"- GEO phenotype rows on Hub pack-membership GSM: **{n_hub_overlap}** (must be 0)",
        f"- Atlas blobs on `sample.metadata_json`: **{n_atlas_on_sample}** (must be 0)",
        f"- `sample_source_membership` includes `geo_metadata_backfill`: "
        f"**{summary['invariants']['sample_source_membership_includes_geo']}** (must be false)",
        "",
        "## Tissue ontology",
        "",
    ]
    tissue_map = stats.get("tissue_map") or {}
    lines.append(
        f"- mapped / unmapped / ambiguous / empty: "
        f"**{tissue_map.get('mapped', 0)}** / **{tissue_map.get('unmapped', 0)}** / "
        f"**{tissue_map.get('ambiguous', 0)}** / **{tissue_map.get('empty', 0)}**"
    )
    unmapped = stats.get("unmapped_tissue_examples") or []
    if unmapped:
        lines.append(f"- Unmapped examples: {', '.join(f'`{t}`' for t in unmapped[:20])}")
    lines.extend(
        [
            "",
            "## Phenotypes by id",
            "",
            "| phenotype_id | rows | observed | unique GSM |",
            "| --- | ---: | ---: | ---: |",
        ]
    )
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
            f"- Unique GSM before (disk census): **{unique_before}**",
            f"- Unique GSM after: **{unique_after}**",
            f"- GEO phenotype rows before merge (in-memory): **{n_geo_before_merge}**",
            f"- GEO observed rows in previous disk census: **{n_geo_obs_before_disk}**",
            f"- GEO observed phenotype rows after: **{geo_after}**",
            f"- Authoritative phenotype rows added this merge: "
            f"**{stats.get('n_phenotype_rows_added', 0)}**",
            "",
        ]
    )
    if dirty_disk_baseline:
        lines.extend(
            [
                "> **Dirty disk baseline:** previous `census.json` already had GEO rows. "
                "That does **not** mean this merge added zero — use "
                "`n_phenotype_rows_added` / in-memory before_merge (0 on a full rebuild).",
                "",
            ]
        )
    lines.extend(
        [
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
