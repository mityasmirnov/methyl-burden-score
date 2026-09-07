"""EWAS DataHub repository census + Hub/EWAS_db lane flags."""

from __future__ import annotations

import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

import pandas as pd

from mbs.annotation.manifest import sha256_file, utc_now_iso
from mbs.platform_id import normalize_platform

DATAHUB_REPO_BASIC = "https://ngdc.cncb.ac.cn/ewas/datahub/repository/basic"
DATAHUB_PLATFORMS = ("450K", "850K", "935K")
OFFICIAL_SAMPLE_TOTAL = 180_317
OFFICIAL_FIELD_TOTAL = 296
OFFICIAL_TISSUE_TOTAL = 1_413
OFFICIAL_DISEASE_TOTAL = 842

CENSUS_PARQUET_NAME = "ewas_datahub_sample_census.parquet"
CENSUS_MANIFEST_NAME = "ewas_datahub_sample_census.manifest.json"
FIELD_VOCAB_NAME = "ewas_datahub_field_vocabulary.json"
TISSUE_VOCAB_NAME = "ewas_datahub_tissue_vocab.parquet"
DISEASE_VOCAB_NAME = "ewas_datahub_disease_vocab.parquet"

CORE_FIELDS = frozenset(
    {
        "sample id",
        "project id",
        "platform",
        "tissue",
        "sample type",
        "disease",
        "age (year)",
        "sex",
        "gender",
    }
)

# DataHub sample type → catalog case_control (never invent control for blanks).
SAMPLE_TYPE_CASE_CONTROL = {
    "control": "control",
    "disease tissue": "case",
    "adjacent normal": "adjacent_normal",
}


def datahub_census_enabled() -> bool:
    return os.environ.get("MBS_SKIP_DATAHUB_CENSUS", "0") != "1"


def catalog_platform_from_datahub(raw: object | None) -> str | None:
    """Map DataHub platform strings (450K/850K/935K) to catalog platform_id."""
    if raw is None or (isinstance(raw, float) and pd.isna(raw)):
        return None
    text = str(raw).strip()
    if not text:
        return None
    aliases = {
        "450K": "HM450",
        "450k": "HM450",
        "850K": "EPIC",
        "850k": "EPIC",
        "935K": "EPICv2",
        "935k": "EPICv2",
        "EPIC": "EPIC",
        "EPICv2": "EPICv2",
        "HM450": "HM450",
    }
    if text in aliases:
        return aliases[text]
    return normalize_platform(text)


def build_sample_lane_flags(
    *,
    samples: pd.DataFrame,
    membership: pd.DataFrame,
    ewas_files: pd.DataFrame,
) -> pd.DataFrame:
    """One row per catalog sample: Hub baseline vs EWAS_db assay membership."""
    if samples.empty:
        return pd.DataFrame(
            columns=[
                "sample_id",
                "study_id",
                "in_hub_baseline",
                "in_ewas_db",
                "hub_families",
            ]
        )
    hub_ids: set[str] = set()
    families_by_sample: dict[str, list[str]] = {}
    if not membership.empty and "sample_id" in membership.columns:
        hub_ids = set(membership["sample_id"].astype(str))
        if "phenotype_family" in membership.columns:
            for sid, grp in membership.groupby(membership["sample_id"].astype(str)):
                fams = sorted(
                    {
                        str(f).strip()
                        for f in grp["phenotype_family"].dropna().astype(str)
                        if str(f).strip()
                    }
                )
                families_by_sample[str(sid)] = fams

    ewas_ids: set[str] = set()
    if not ewas_files.empty and "sample_id" in ewas_files.columns:
        ewas_ids = set(ewas_files["sample_id"].astype(str))

    rows: list[dict[str, Any]] = []
    for rec in samples.to_dict(orient="records"):
        sid = str(rec["sample_id"])
        in_hub = sid in hub_ids
        in_ewas = sid in ewas_ids
        rows.append(
            {
                "sample_id": sid,
                "study_id": str(rec.get("study_id") or ""),
                "in_hub_baseline": bool(in_hub),
                "in_ewas_db": bool(in_ewas),
                "hub_families": json.dumps(families_by_sample.get(sid, [])),
            }
        )
    return pd.DataFrame(rows)


def stamp_lane_flags_into_metadata(
    *,
    samples: pd.DataFrame,
    studies: pd.DataFrame,
    sample_flags: pd.DataFrame,
    study_flags: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Mirror lane flags into metadata_json.lanes (alongside first-class tables)."""
    sample_by = (
        sample_flags.set_index("sample_id", drop=False)
        if not sample_flags.empty
        else pd.DataFrame()
    )
    study_by = (
        study_flags.set_index("study_id", drop=False)
        if not study_flags.empty
        else pd.DataFrame()
    )
    out_samples: list[dict[str, Any]] = []
    for rec in samples.to_dict(orient="records"):
        sid = str(rec["sample_id"])
        meta: dict[str, Any] = {}
        raw = rec.get("metadata_json")
        if raw:
            meta = json.loads(raw) if isinstance(raw, str) else dict(raw)
        lanes: list[str] = []
        if sid in sample_by.index:
            row = sample_by.loc[sid]
            if bool(row["in_hub_baseline"]):
                lanes.append("ewas_datahub_baseline")
            if bool(row["in_ewas_db"]):
                lanes.append("ewas_datahub_db")
        meta["lanes"] = lanes
        rec["metadata_json"] = json.dumps(meta, sort_keys=True)
        out_samples.append(rec)

    out_studies: list[dict[str, Any]] = []
    for rec in studies.to_dict(orient="records"):
        study_id = str(rec["study_id"])
        meta: dict[str, Any] = {}
        raw = rec.get("metadata_json")
        if raw:
            meta = json.loads(raw) if isinstance(raw, str) else dict(raw)
        lanes: list[str] = []
        if study_id in study_by.index:
            row = study_by.loc[study_id]
            if bool(row["in_hub_baseline"]):
                lanes.append("ewas_datahub_baseline")
            if bool(row["in_ewas_db"]):
                lanes.append("ewas_datahub_db")
        # Preserve any other pre-existing lane tags.
        existing = meta.get("lanes") or []
        if isinstance(existing, list):
            for lane in existing:
                if lane not in lanes:
                    lanes.append(str(lane))
        meta["lanes"] = lanes
        rec["metadata_json"] = json.dumps(meta, sort_keys=True)
        out_studies.append(rec)

    return pd.DataFrame(out_samples), pd.DataFrame(out_studies)


def build_study_lane_flags(
    *,
    studies: pd.DataFrame,
    sample_flags: pd.DataFrame,
) -> pd.DataFrame:
    """Aggregate sample lane flags to study level."""
    if studies.empty:
        return pd.DataFrame(
            columns=[
                "study_id",
                "in_hub_baseline",
                "in_ewas_db",
                "n_hub_samples",
                "n_ewas_db_samples",
                "n_samples",
            ]
        )
    if sample_flags.empty:
        return pd.DataFrame(
            [
                {
                    "study_id": str(rec["study_id"]),
                    "in_hub_baseline": False,
                    "in_ewas_db": False,
                    "n_hub_samples": 0,
                    "n_ewas_db_samples": 0,
                    "n_samples": 0,
                }
                for rec in studies.to_dict(orient="records")
            ]
        )
    agg = (
        sample_flags.groupby("study_id", as_index=False)
        .agg(
            n_samples=("sample_id", "count"),
            n_hub_samples=("in_hub_baseline", "sum"),
            n_ewas_db_samples=("in_ewas_db", "sum"),
        )
        .astype({"n_hub_samples": "int64", "n_ewas_db_samples": "int64", "n_samples": "int64"})
    )
    agg["in_hub_baseline"] = agg["n_hub_samples"] > 0
    agg["in_ewas_db"] = agg["n_ewas_db_samples"] > 0
    by_study = agg.set_index("study_id", drop=False)
    rows: list[dict[str, Any]] = []
    for rec in studies.to_dict(orient="records"):
        sid = str(rec["study_id"])
        if sid in by_study.index:
            row = by_study.loc[sid]
            rows.append(
                {
                    "study_id": sid,
                    "in_hub_baseline": bool(row["in_hub_baseline"]),
                    "in_ewas_db": bool(row["in_ewas_db"]),
                    "n_hub_samples": int(row["n_hub_samples"]),
                    "n_ewas_db_samples": int(row["n_ewas_db_samples"]),
                    "n_samples": int(row["n_samples"]),
                }
            )
        else:
            rows.append(
                {
                    "study_id": sid,
                    "in_hub_baseline": False,
                    "in_ewas_db": False,
                    "n_hub_samples": 0,
                    "n_ewas_db_samples": 0,
                    "n_samples": 0,
                }
            )
    return pd.DataFrame(rows)


def _http_get_json(url: str, *, timeout_s: float = 120.0) -> dict[str, Any]:
    if not url.startswith(("http://", "https://")):
        raise ValueError(f"refusing non-http(s) URL: {url!r}")
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": "methyl-burden-score/1.0 (research; datahub census)",
            "Accept": "application/json, text/javascript, */*; q=0.01",
            "X-Requested-With": "XMLHttpRequest",
            "Referer": "https://ngdc.cncb.ac.cn/ewas/datahub/repository",
        },
    )
    last_exc: BaseException | None = None
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=timeout_s) as resp:  # noqa: S310
                return json.loads(resp.read().decode("utf-8"))
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            last_exc = exc
            time.sleep(min(30.0, 1.5 * (2**attempt)))
    assert last_exc is not None
    raise last_exc


def fetch_platform_page(
    *,
    platform: str,
    offset: int,
    limit: int,
) -> dict[str, Any]:
    params = urllib.parse.urlencode(
        {
            "field": "platform",
            "val": platform,
            "offset": str(offset),
            "limit": str(limit),
            "relationship": "",
        },
        doseq=True,
    )
    return _http_get_json(f"{DATAHUB_REPO_BASIC}?{params}")


def fetch_datahub_repository_census(
    *,
    cache_root: Path,
    output_dir: Path,
    page_size: int = 200,
    delay_s: float = 0.2,
    platforms: tuple[str, ...] = DATAHUB_PLATFORMS,
    from_cache_only: bool = False,
) -> dict[str, Any]:
    """Paginate repository/basic by platform; write census Parquet + vocabs."""
    cache_root = cache_root.resolve()
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    all_rows: list[dict[str, Any]] = []
    platform_totals: dict[str, int] = {}
    field_keys: set[str] = set()

    for platform in platforms:
        plat_cache = cache_root / platform
        plat_cache.mkdir(parents=True, exist_ok=True)
        offset = 0
        total: int | None = None
        while True:
            page_path = plat_cache / f"page_{offset:06d}.json"
            if page_path.is_file():
                payload = json.loads(page_path.read_text(encoding="utf-8"))
            elif from_cache_only:
                break
            else:
                try:
                    payload = fetch_platform_page(
                        platform=platform, offset=offset, limit=page_size
                    )
                except (urllib.error.URLError, TimeoutError, OSError, json.JSONDecodeError) as exc:
                    raise RuntimeError(
                        f"DataHub repository fetch failed "
                        f"platform={platform} offset={offset}: {exc}"
                    ) from exc
                page_path.write_text(json.dumps(payload), encoding="utf-8")
                time.sleep(delay_s)
            if total is None:
                total = int(str(payload.get("total") or 0))
                platform_totals[platform] = total
            content = payload.get("content") or []
            if not content:
                break
            for row in content:
                if not isinstance(row, dict):
                    continue
                cleaned = {str(k): v for k, v in row.items() if k != "0"}
                field_keys.update(cleaned.keys())
                cleaned["_platform_query"] = platform
                all_rows.append(cleaned)
            offset += page_size
            if total is not None and offset >= total:
                break
            if len(content) < page_size:
                break

    if not all_rows:
        raise RuntimeError("DataHub census fetch produced zero rows")

    frame = pd.DataFrame(all_rows)
    # Normalize core columns
    rename = {
        "sample id": "sample_id",
        "project id": "project_id",
        "sample type": "sample_type",
        "age (year)": "age_year",
    }
    for old, new in rename.items():
        if old in frame.columns and new not in frame.columns:
            frame[new] = frame[old]
    if "sample_id" not in frame.columns:
        raise RuntimeError("census rows missing sample id")
    frame["sample_id"] = frame["sample_id"].astype(str)
    frame = frame.drop_duplicates(subset=["sample_id"], keep="first")
    if "platform" in frame.columns:
        frame["catalog_platform_id"] = [
            catalog_platform_from_datahub(v) for v in frame["platform"].tolist()
        ]
    else:
        frame["catalog_platform_id"] = None

    census_path = output_dir / CENSUS_PARQUET_NAME
    frame.to_parquet(census_path, index=False)

    # Field vocabulary
    field_vocab = {
        "n_fields": len(field_keys),
        "official_n_fields": OFFICIAL_FIELD_TOTAL,
        "fields": sorted(field_keys),
        "core_fields_present": sorted(field_keys & CORE_FIELDS),
    }
    (output_dir / FIELD_VOCAB_NAME).write_text(
        json.dumps(field_vocab, indent=2) + "\n", encoding="utf-8"
    )

    tissue_vocab = _value_counts_frame(frame, "tissue")
    disease_vocab = _value_counts_frame(frame, "disease")
    tissue_vocab.to_parquet(output_dir / TISSUE_VOCAB_NAME, index=False)
    disease_vocab.to_parquet(output_dir / DISEASE_VOCAB_NAME, index=False)

    manifest = {
        "generated_at": utc_now_iso(),
        "n_samples": len(frame),
        "official_n_samples": OFFICIAL_SAMPLE_TOTAL,
        "platform_totals_api": platform_totals,
        "platform_counts_local": (
            {str(k): int(v) for k, v in frame["platform"].astype(str).value_counts().items()}
            if "platform" in frame.columns
            else {}
        ),
        "n_tissues": len(tissue_vocab),
        "n_diseases": len(disease_vocab),
        "official_n_tissues": OFFICIAL_TISSUE_TOTAL,
        "official_n_diseases": OFFICIAL_DISEASE_TOTAL,
        "n_fields": field_vocab["n_fields"],
        "official_n_fields": OFFICIAL_FIELD_TOTAL,
        "census_parquet": str(census_path),
        "census_sha256": sha256_file(census_path),
        "cache_root": str(cache_root),
    }
    (output_dir / CENSUS_MANIFEST_NAME).write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def _value_counts_frame(frame: pd.DataFrame, column: str) -> pd.DataFrame:
    if column not in frame.columns:
        return pd.DataFrame(columns=["value", "n"])
    return (
        frame[column]
        .dropna()
        .astype(str)
        .str.strip()
        .loc[lambda s: s != ""]
        .value_counts()
        .rename_axis("value")
        .reset_index(name="n")
    )

def load_datahub_census(data_root: Path) -> pd.DataFrame:
    path = data_root / "canonical" / "phenotypes" / CENSUS_PARQUET_NAME
    if not path.is_file():
        return pd.DataFrame()
    return pd.read_parquet(path)


def merge_datahub_census(
    *,
    samples: pd.DataFrame,
    phenotypes: pd.DataFrame,
    studies: pd.DataFrame,
    census: pd.DataFrame,
    hub_sample_ids: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Merge DataHub census into catalog tables.

    Hub pack samples keep existing phenotypes; census fills platform/tissue/
    sample_type gaps and attaches a compact ``metadata_json.datahub`` bag for
    matched GSMs. Full sparse field dumps stay in the census Parquet (not
    copied into every sample JSON).
    """
    stats: dict[str, Any] = {
        "n_census_rows": len(census),
        "n_samples_matched": 0,
        "n_platform_filled": 0,
        "n_tissue_filled": 0,
        "n_sample_type_filled": 0,
        "n_phenotype_rows_added": 0,
        "n_studies_platform_set": 0,
    }
    if samples.empty or census.empty or "sample_id" not in census.columns:
        return samples, phenotypes, studies, stats

    census_u = census.drop_duplicates(subset=["sample_id"]).copy()
    census_u["sample_id"] = census_u["sample_id"].astype(str)
    if "sample_type" not in census_u.columns and "sample type" in census_u.columns:
        census_u["sample_type"] = census_u["sample type"]
    if "age_year" not in census_u.columns and "age (year)" in census_u.columns:
        census_u["age_year"] = census_u["age (year)"]
    if "catalog_platform_id" not in census_u.columns:
        census_u["catalog_platform_id"] = [
            catalog_platform_from_datahub(v) for v in census_u.get("platform", pd.Series(dtype=object))
        ]
    else:
        missing_plat = census_u["catalog_platform_id"].map(_blank)
        if missing_plat.any() and "platform" in census_u.columns:
            census_u.loc[missing_plat, "catalog_platform_id"] = [
                catalog_platform_from_datahub(v)
                for v in census_u.loc[missing_plat, "platform"].tolist()
            ]

    samples_out = samples.copy()
    samples_out["sample_id"] = samples_out["sample_id"].astype(str)
    joined = samples_out.merge(
        census_u,
        on="sample_id",
        how="left",
        suffixes=("", "_census"),
    )
    matched_mask = joined["platform"].notna() if "platform" in joined.columns else joined["catalog_platform_id"].notna()
    # If platform column missing, use catalog_platform_id / project_id
    if "platform" not in joined.columns:
        matched_mask = joined["catalog_platform_id"].notna()
    else:
        # Rows present in census: merge keeps census columns; unmatched have NA platform
        matched_mask = joined["sample_id"].isin(set(census_u["sample_id"]))
    stats["n_samples_matched"] = int(matched_mask.sum())

    tissue_blank = joined["tissue_raw"].map(_blank) if "tissue_raw" in joined.columns else pd.Series(True, index=joined.index)
    if "tissue" in joined.columns:
        fill_tissue = matched_mask & tissue_blank & ~joined["tissue"].map(_blank)
        joined.loc[fill_tissue, "tissue_raw"] = joined.loc[fill_tissue, "tissue"].astype(str).str.strip()
        stats["n_tissue_filled"] = int(fill_tissue.sum())

    cc_blank = joined["case_control"].map(_blank) if "case_control" in joined.columns else pd.Series(True, index=joined.index)
    if "sample_type" in joined.columns:
        mapped_cc = joined["sample_type"].map(
            lambda v: SAMPLE_TYPE_CASE_CONTROL.get(str(v).strip().lower()) if not _blank(v) else None
        )
        fill_cc = matched_mask & cc_blank & mapped_cc.notna()
        joined.loc[fill_cc, "case_control"] = mapped_cc.loc[fill_cc]
        stats["n_sample_type_filled"] = int(fill_cc.sum())

    plat_ok = matched_mask & ~joined["catalog_platform_id"].map(_blank)
    stats["n_platform_filled"] = int(plat_ok.sum())

    # Compact datahub bag on matched rows only (vectorized apply on subset).
    bag_cols = [
        c
        for c in (
            "platform",
            "tissue",
            "sample_type",
            "disease",
            "age_year",
            "sex",
            "gender",
            "project_id",
            "catalog_platform_id",
        )
        if c in joined.columns
    ]

    def _merge_meta(row: pd.Series) -> str | None:
        if not bool(row.get("_matched")):
            raw = row.get("metadata_json")
            return None if (raw is None or (isinstance(raw, float) and pd.isna(raw))) else str(raw)
        meta: dict[str, Any] = {}
        raw = row.get("metadata_json")
        if raw and not (isinstance(raw, float) and pd.isna(raw)):
            meta = json.loads(raw) if isinstance(raw, str) else dict(raw)
        blob: dict[str, Any] = {}
        for col in bag_cols:
            val = row.get(col)
            if _blank(val):
                continue
            if col in {"age_year", "catalog_platform_id"} or isinstance(val, (int, float)):
                if col == "age_year":
                    parsed = _parse_float(val)
                    if parsed is not None:
                        blob[col] = parsed
                else:
                    blob[col] = str(val).strip() if col == "catalog_platform_id" else val
            else:
                blob[col] = str(val).strip()
        if blob:
            meta["datahub"] = blob
        return json.dumps(meta, sort_keys=True) if meta else (
            None if raw is None or (isinstance(raw, float) and pd.isna(raw)) else str(raw)
        )

    joined["_matched"] = matched_mask
    if matched_mask.any():
        joined.loc[matched_mask, "metadata_json"] = joined.loc[matched_mask].apply(
            _merge_meta, axis=1
        )

    # Drop census-only columns before returning samples schema.
    sample_cols = list(samples.columns)
    samples_out = joined[sample_cols].copy()

    # Phenotype rows for non-Hub matched samples.
    non_hub = matched_mask & ~joined["sample_id"].isin(hub_sample_ids)
    new_pheno_rows: list[dict[str, Any]] = []
    source_family = "ewas_datahub_repository"
    if non_hub.any():
        sub = joined.loc[non_hub]
        if "tissue" in sub.columns:
            for sid, tissue in zip(
                sub["sample_id"], sub["tissue"], strict=True
            ):
                if not _blank(tissue):
                    new_pheno_rows.append(
                        _pheno_row(
                            str(sid),
                            "tissue",
                            categorical=str(tissue).strip(),
                            source_family=source_family,
                        )
                    )
        sex_series = sub["sex"] if "sex" in sub.columns else None
        gender_series = sub["gender"] if "gender" in sub.columns else None
        for i, sid in enumerate(sub["sample_id"].tolist()):
            sex_val = None
            if sex_series is not None:
                sex_val = sex_series.iloc[i]
            if _blank(sex_val) and gender_series is not None:
                sex_val = gender_series.iloc[i]
            if not _blank(sex_val):
                new_pheno_rows.append(
                    _pheno_row(
                        str(sid),
                        "sex",
                        categorical=str(sex_val).strip().lower(),
                        source_family=source_family,
                    )
                )
        if "age_year" in sub.columns:
            for sid, age_val in zip(sub["sample_id"], sub["age_year"], strict=True):
                age_num = _parse_float(age_val)
                if age_num is not None and 0 <= age_num <= 120:
                    new_pheno_rows.append(
                        _pheno_row(str(sid), "age", numeric=age_num, source_family=source_family)
                    )
        if "disease" in sub.columns and "sample_type" in sub.columns:
            for sid, disease, sample_type in zip(
                sub["sample_id"], sub["disease"], sub["sample_type"], strict=True
            ):
                if not _blank(disease) and mapped_case_control_is_case(sample_type):
                    new_pheno_rows.append(
                        _pheno_row(
                            str(sid),
                            "disease",
                            categorical=str(disease).strip(),
                            source_family=source_family,
                            label_status="case",
                        )
                    )

    if new_pheno_rows:
        add = pd.DataFrame(new_pheno_rows)
        stats["n_phenotype_rows_added"] = len(add)
        phenotypes = (
            pd.concat([phenotypes, add], ignore_index=True)
            if not phenotypes.empty
            else add
        )
        phenotypes = phenotypes.drop_duplicates(
            subset=["sample_id", "phenotype_id", "source_family"], keep="first"
        )

    # Study platform_id when a single catalog platform covers the study's matched samples.
    if not studies.empty:
        matched_samples = joined.loc[matched_mask, ["study_id", "catalog_platform_id"]].copy()
        matched_samples["study_id"] = matched_samples["study_id"].astype(str)
        matched_samples = matched_samples[~matched_samples["catalog_platform_id"].map(_blank)]
        study_plats: dict[str, set[str]] = {}
        if not matched_samples.empty:
            for study_id, grp in matched_samples.groupby("study_id"):
                plats = {str(p).strip() for p in grp["catalog_platform_id"].tolist() if str(p).strip()}
                if plats:
                    study_plats[str(study_id)] = plats
        study_rows = []
        for rec in studies.to_dict(orient="records"):
            study_id = str(rec["study_id"])
            plats = study_plats.get(study_id) or set()
            if _blank(rec.get("platform_id")) and len(plats) == 1:
                rec["platform_id"] = next(iter(plats))
                stats["n_studies_platform_set"] += 1
            meta: dict[str, Any] = {}
            raw = rec.get("metadata_json")
            if raw:
                meta = json.loads(raw) if isinstance(raw, str) else dict(raw)
            if plats:
                meta.setdefault("datahub", {})
                if isinstance(meta["datahub"], dict):
                    meta["datahub"]["platforms"] = sorted(plats)
                rec["metadata_json"] = json.dumps(meta, sort_keys=True)
            study_rows.append(rec)
        studies = pd.DataFrame(study_rows)

    return samples_out, phenotypes, studies, stats


def mapped_case_control_is_case(sample_type: object) -> bool:
    if _blank(sample_type):
        return False
    return SAMPLE_TYPE_CASE_CONTROL.get(str(sample_type).strip().lower()) == "case"


def _pheno_row(
    sample_id: str,
    phenotype_id: str,
    *,
    source_family: str,
    categorical: str | None = None,
    numeric: float | None = None,
    label_status: str = "observed",
) -> dict[str, Any]:
    return {
        "sample_id": sample_id,
        "phenotype_id": phenotype_id,
        "numeric_value": numeric,
        "categorical_value": categorical,
        "label_status": label_status,
        "is_observed": True,
        "source_family": source_family,
        "source_record_id": None,
        "ontology_id": None,
    }


def _blank(value: object) -> bool:
    if value is None:
        return True
    if isinstance(value, float) and pd.isna(value):
        return True
    return str(value).strip() == "" or str(value).strip().lower() == "nan"


def _parse_float(value: object) -> float | None:
    if _blank(value):
        return None
    try:
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None


def write_lane_flags_report(
    *,
    sample_flags: pd.DataFrame,
    study_flags: pd.DataFrame,
    report_dir: Path,
) -> Path:
    """Write Hub vs EWAS_db membership census."""
    report_dir.mkdir(parents=True, exist_ok=True)
    n = len(sample_flags)
    in_hub = int(sample_flags["in_hub_baseline"].sum()) if n else 0
    in_ewas = int(sample_flags["in_ewas_db"].sum()) if n else 0
    both = (
        int(((sample_flags["in_hub_baseline"]) & (sample_flags["in_ewas_db"])).sum()) if n else 0
    )
    hub_only = in_hub - both
    ewas_only = in_ewas - both
    neither = n - hub_only - ewas_only - both

    sn = len(study_flags)
    s_hub = int(study_flags["in_hub_baseline"].sum()) if sn else 0
    s_ewas = int(study_flags["in_ewas_db"].sum()) if sn else 0
    s_both = (
        int(((study_flags["in_hub_baseline"]) & (study_flags["in_ewas_db"])).sum()) if sn else 0
    )

    payload = {
        "n_samples": n,
        "samples": {
            "hub_only": hub_only,
            "ewas_db_only": ewas_only,
            "both": both,
            "neither": neither,
            "in_hub_baseline": in_hub,
            "in_ewas_db": in_ewas,
        },
        "n_studies": sn,
        "studies": {
            "hub_only": s_hub - s_both,
            "ewas_db_only": s_ewas - s_both,
            "both": s_both,
            "in_hub_baseline": s_hub,
            "in_ewas_db": s_ewas,
        },
        "note": (
            "Hub baseline = sample_source_membership (nine packs). "
            "EWAS_db = assay .txt on disk. Training phenotypes still come from Hub packs."
        ),
    }
    json_path = report_dir / "hub_vs_ewas_db_membership.json"
    json_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    lines = [
        "# Hub vs EWAS_db membership",
        "",
        f"- Samples: **{n}**",
        f"  - Hub-only: **{hub_only}**",
        f"  - EWAS_db-only: **{ewas_only}**",
        f"  - Both: **{both}**",
        f"  - Neither: **{neither}**",
        "",
        f"- Studies: **{sn}**",
        f"  - Hub-only: **{s_hub - s_both}**",
        f"  - EWAS_db-only: **{s_ewas - s_both}**",
        f"  - Both: **{s_both}**",
        "",
        payload["note"],
        "",
    ]
    (report_dir / "hub_vs_ewas_db_membership.md").write_text("\n".join(lines), encoding="utf-8")
    return json_path


def write_datahub_census_report(
    *,
    census: pd.DataFrame,
    sample_flags: pd.DataFrame,
    merge_stats: dict[str, Any],
    report_dir: Path,
    manifest: dict[str, Any] | None = None,
) -> Path:
    """Compare official DataHub totals to local census + catalog coverage."""
    out = report_dir / "datahub_metadata_census"
    out.mkdir(parents=True, exist_ok=True)
    n_local = len(census) if not census.empty else 0
    plat_counts = (
        {str(k): int(v) for k, v in census["platform"].astype(str).value_counts().items()}
        if not census.empty and "platform" in census.columns
        else {}
    )
    n_tissues = (
        int(census["tissue"].dropna().astype(str).str.strip().loc[lambda s: s != ""].nunique())
        if not census.empty and "tissue" in census.columns
        else 0
    )
    n_diseases = (
        int(census["disease"].dropna().astype(str).str.strip().loc[lambda s: s != ""].nunique())
        if not census.empty and "disease" in census.columns
        else 0
    )
    catalog_ids = set(sample_flags["sample_id"].astype(str)) if not sample_flags.empty else set()
    census_ids = set(census["sample_id"].astype(str)) if not census.empty else set()
    in_both = len(catalog_ids & census_ids)
    census_not_catalog = len(census_ids - catalog_ids)
    catalog_not_census = len(catalog_ids - census_ids)
    hub_in_census = 0
    if not sample_flags.empty and census_ids:
        hub_mask = sample_flags["in_hub_baseline"].astype(bool)
        hub_in_census = int(
            sample_flags.loc[hub_mask, "sample_id"].astype(str).isin(census_ids).sum()
        )

    payload = {
        "official": {
            "n_samples": OFFICIAL_SAMPLE_TOTAL,
            "n_fields": OFFICIAL_FIELD_TOTAL,
            "n_tissues": OFFICIAL_TISSUE_TOTAL,
            "n_diseases": OFFICIAL_DISEASE_TOTAL,
        },
        "local_census": {
            "n_samples": n_local,
            "platform_counts": plat_counts,
            "n_tissues": n_tissues,
            "n_diseases": n_diseases,
            "n_fields": (manifest or {}).get("n_fields"),
        },
        "vs_catalog": {
            "n_overlap": in_both,
            "census_not_in_catalog": census_not_catalog,
            "catalog_not_in_census": catalog_not_census,
            "hub_samples_in_census": hub_in_census,
        },
        "merge_stats": merge_stats,
        "manifest": manifest,
    }
    json_path = out / "summary.json"
    json_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    lines = [
        "# DataHub metadata census",
        "",
        "## Official vs local census Parquet",
        "",
        "| Metric | Official | Local |",
        "| --- | ---: | ---: |",
        f"| Samples | {OFFICIAL_SAMPLE_TOTAL} | {n_local} |",
        f"| Tissues/Cells | {OFFICIAL_TISSUE_TOTAL} | {n_tissues} |",
        f"| Diseases | {OFFICIAL_DISEASE_TOTAL} | {n_diseases} |",
        f"| Fields | {OFFICIAL_FIELD_TOTAL} | {(manifest or {}).get('n_fields', '—')} |",
        "",
        "## Platform split (local census)",
        "",
    ]
    for plat, count in sorted(plat_counts.items()):
        lines.append(f"- `{plat}`: {count}")
    lines.extend(
        [
            "",
            "## Catalog coverage",
            "",
            f"- Overlap census ∩ catalog: **{in_both}**",
            f"- In census, not yet in catalog (missing betas): **{census_not_catalog}**",
            f"- In catalog, not in census: **{catalog_not_census}**",
            f"- Hub baseline samples present in census: **{hub_in_census}**",
            "",
            "## Merge into refresh",
            "",
            f"- Enabled stats: `{json.dumps(merge_stats)}`",
            "",
            "Hub packs remain training phenotype SoT; DataHub census fills platform/"
            "tissue/sample_type and `metadata_json.datahub` for stratification.",
            "",
        ]
    )
    (out / "summary.md").write_text("\n".join(lines), encoding="utf-8")
    return json_path
