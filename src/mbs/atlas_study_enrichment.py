"""Study-level EWAS Atlas enrichment for catalog stratification (not sample labels)."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pandas as pd

from mbs.inspect_ewas_metadata import read_atlas_tsv

DEFAULT_ATLAS_GSE_MAP = Path("configs/data/atlas_gse_es_map.tsv")
ATLAS_STUDY_ID_RE = re.compile(r"^ES\d+$", re.IGNORECASE)


def gds_uid_for_gse(gse_id: str) -> str:
    """NCBI GEO DataSets uid for a series accession (200 + zero-padded GSE number)."""
    suffix = gse_id.strip().upper().removeprefix("GSE")
    if not suffix.isdigit():
        raise ValueError(f"not a GSE accession: {gse_id!r}")
    return f"200{int(suffix):06d}"


def _sorted_unique(values: pd.Series) -> list[str]:
    cleaned = [str(v).strip() for v in values.dropna().astype(str) if str(v).strip()]
    return sorted(dict.fromkeys(cleaned))


def _parse_int(value: object) -> int | None:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def load_atlas_gse_map(path: Path) -> pd.DataFrame:
    """Optional curated GSE↔Atlas map (never join GSE to ES* by string equality)."""
    columns = ["gse_id", "atlas_study_id", "pmid", "source"]
    if not path.is_file():
        return pd.DataFrame(columns=columns)
    try:
        frame = pd.read_csv(path, sep="\t", dtype=str, comment="#")
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=columns)
    if frame.empty:
        return pd.DataFrame(columns=columns)
    frame = frame.fillna("")
    for col in columns:
        if col not in frame.columns:
            frame[col] = ""
    frame = frame[columns].copy()
    frame["gse_id"] = frame["gse_id"].str.strip()
    frame["atlas_study_id"] = frame["atlas_study_id"].str.strip()
    frame["pmid"] = frame["pmid"].str.strip()
    frame["source"] = frame["source"].str.strip()
    return frame.loc[frame["gse_id"] != ""].reset_index(drop=True)


def build_atlas_reference(
    *,
    atlas_root: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, list[str]]]:
    """Load Atlas studies/cohorts and PMID→atlas_study_id index."""
    studies_path = atlas_root / "EWAS_Atlas_studies.tsv"
    cohorts_path = atlas_root / "EWAS_Atlas_cohorts.tsv"
    if not studies_path.is_file() or not cohorts_path.is_file():
        return pd.DataFrame(), pd.DataFrame(), {}
    studies = read_atlas_tsv(studies_path)
    cohorts = read_atlas_tsv(cohorts_path)
    pmid_to_study: dict[str, list[str]] = {}
    if not studies.empty and {"study_ID", "PMID"}.issubset(studies.columns):
        for _, row in studies.drop_duplicates(subset=["study_ID", "PMID"]).iterrows():
            pmid = str(row["PMID"]).strip()
            study_id = str(row["study_ID"]).strip()
            if not pmid or not study_id:
                continue
            pmid_to_study.setdefault(pmid, [])
            if study_id not in pmid_to_study[pmid]:
                pmid_to_study[pmid].append(study_id)
        for pmid in pmid_to_study:
            pmid_to_study[pmid] = sorted(pmid_to_study[pmid])
    return studies, cohorts, pmid_to_study


def _summarize_atlas_studies(
    *,
    atlas_study_ids: list[str],
    studies: pd.DataFrame,
    cohorts: pd.DataFrame,
) -> dict[str, Any]:
    atlas_set = set(atlas_study_ids)
    study_slice = studies.loc[studies["study_ID"].astype(str).isin(atlas_set)]
    cohort_slice = cohorts.loc[cohorts["study_ID"].astype(str).isin(atlas_set)]

    traits = _sorted_unique(study_slice["trait"]) if "trait" in study_slice.columns else []
    pmids = _sorted_unique(study_slice["PMID"]) if "PMID" in study_slice.columns else []

    cohort_descriptions: list[str] = []
    if "description" in cohort_slice.columns:
        cohort_descriptions = _sorted_unique(cohort_slice["description"])
    tissues = _sorted_unique(cohort_slice["tissue"]) if "tissue" in cohort_slice.columns else []
    platforms = _sorted_unique(cohort_slice["platform"]) if "platform" in cohort_slice.columns else []
    ancestries = (
        _sorted_unique(cohort_slice["ancestry"]) if "ancestry" in cohort_slice.columns else []
    )

    sample_sizes = [
        n
        for n in (_parse_int(v) for v in cohort_slice.get("sample_size", pd.Series(dtype=object)))
        if n is not None
    ]
    return {
        "atlas_study_ids": sorted(atlas_set),
        "pmid": pmids[0] if len(pmids) == 1 else None,
        "pmids": pmids,
        "atlas_traits": traits,
        "n_atlas_cohorts": int(len(cohort_slice)),
        "total_sample_size": int(sum(sample_sizes)) if sample_sizes else None,
        "tissues": tissues,
        "cohort_descriptions": cohort_descriptions,
        "platforms": platforms,
        "ancestries": ancestries,
    }


def resolve_atlas_study_ids(
    *,
    catalog_study_id: str,
    gse_map: pd.DataFrame,
    pmid_to_study: dict[str, list[str]],
) -> tuple[list[str], str, str | None]:
    """Return (atlas_study_ids, join_method, pmid_used)."""
    study_id = catalog_study_id.strip()
    if ATLAS_STUDY_ID_RE.match(study_id):
        return [study_id.upper()], "atlas_study_id", None

    map_row = gse_map.loc[gse_map["gse_id"] == study_id]
    if not map_row.empty:
        row = map_row.iloc[0]
        atlas_ids = [
            x.strip()
            for x in str(row.get("atlas_study_id", "")).split(";")
            if x.strip()
        ]
        pmid = str(row.get("pmid", "")).strip() or None
        if atlas_ids:
            return sorted(dict.fromkeys(atlas_ids)), "gse_es_map", pmid
        if pmid and pmid in pmid_to_study:
            return pmid_to_study[pmid], "pmid", pmid

    return [], "none", None


def build_study_atlas_enrichment(
    *,
    catalog_study_ids: list[str],
    atlas_root: Path,
    gse_map_path: Path,
) -> pd.DataFrame:
    """One row per catalog study with Atlas context for external stratification."""
    studies, cohorts, pmid_to_study = build_atlas_reference(atlas_root=atlas_root)
    if studies.empty or cohorts.empty:
        return pd.DataFrame(
            columns=[
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
            ]
        )

    gse_map = load_atlas_gse_map(gse_map_path)
    rows: list[dict[str, Any]] = []
    for catalog_study_id in sorted(set(catalog_study_ids)):
        atlas_ids, join_method, pmid_used = resolve_atlas_study_ids(
            catalog_study_id=catalog_study_id,
            gse_map=gse_map,
            pmid_to_study=pmid_to_study,
        )
        if not atlas_ids:
            rows.append(
                {
                    "study_id": catalog_study_id,
                    "join_method": join_method,
                    "atlas_study_ids": json.dumps([]),
                    "pmid": None,
                    "n_atlas_cohorts": None,
                    "total_sample_size": None,
                    "tissues": json.dumps([]),
                    "cohort_descriptions": json.dumps([]),
                    "platforms": json.dumps([]),
                    "ancestries": json.dumps([]),
                    "atlas_traits": json.dumps([]),
                }
            )
            continue
        summary = _summarize_atlas_studies(
            atlas_study_ids=atlas_ids,
            studies=studies,
            cohorts=cohorts,
        )
        pmid = pmid_used or summary.get("pmid")
        rows.append(
            {
                "study_id": catalog_study_id,
                "join_method": join_method,
                "atlas_study_ids": json.dumps(summary["atlas_study_ids"]),
                "pmid": pmid,
                "n_atlas_cohorts": summary["n_atlas_cohorts"],
                "total_sample_size": summary["total_sample_size"],
                "tissues": json.dumps(summary["tissues"]),
                "cohort_descriptions": json.dumps(summary["cohort_descriptions"]),
                "platforms": json.dumps(summary["platforms"]),
                "ancestries": json.dumps(summary["ancestries"]),
                "atlas_traits": json.dumps(summary["atlas_traits"]),
            }
        )
    return pd.DataFrame(rows)


def merge_atlas_enrichment_into_studies(
    studies: pd.DataFrame,
    enrichment: pd.DataFrame,
) -> pd.DataFrame:
    """Attach atlas_enrichment blob to study.metadata_json (stratification only)."""
    if studies.empty or enrichment.empty:
        return studies
    by_study = enrichment.set_index("study_id", drop=False)
    updated: list[dict[str, Any]] = []
    for rec in studies.to_dict(orient="records"):
        study_id = str(rec["study_id"])
        meta: dict[str, Any] = {}
        raw = rec.get("metadata_json")
        if raw:
            meta = json.loads(raw) if isinstance(raw, str) else dict(raw)
        if study_id not in by_study.index:
            updated.append(rec)
            continue
        row = by_study.loc[study_id]
        if row["join_method"] == "none":
            updated.append(rec)
            continue
        meta["atlas_enrichment"] = {
            "join_method": str(row["join_method"]),
            "atlas_study_ids": json.loads(str(row["atlas_study_ids"])),
            "pmid": None if pd.isna(row["pmid"]) else str(row["pmid"]),
            "n_atlas_cohorts": None
            if pd.isna(row["n_atlas_cohorts"])
            else int(row["n_atlas_cohorts"]),
            "total_sample_size": None
            if pd.isna(row["total_sample_size"])
            else int(row["total_sample_size"]),
            "tissues": json.loads(str(row["tissues"])),
            "cohort_descriptions": json.loads(str(row["cohort_descriptions"])),
            "platforms": json.loads(str(row["platforms"])),
            "ancestries": json.loads(str(row["ancestries"])),
            "atlas_traits": json.loads(str(row["atlas_traits"])),
        }
        rec = dict(rec)
        rec["metadata_json"] = json.dumps(meta, sort_keys=True)
        updated.append(rec)
    return pd.DataFrame(updated)


def write_study_atlas_enrichment_report(
    *,
    enrichment: pd.DataFrame,
    report_dir: Path,
) -> Path:
    """Write JSON + markdown summary for inspection (no sample labels)."""
    report_dir.mkdir(parents=True, exist_ok=True)
    matched = enrichment.loc[enrichment["join_method"] != "none"] if not enrichment.empty else enrichment
    by_method = (
        matched["join_method"].value_counts().to_dict() if not matched.empty else {}
    )
    payload = {
        "n_catalog_studies": int(len(enrichment)),
        "n_matched": int(len(matched)),
        "by_join_method": by_method,
        "note": (
            "Study-level Atlas context for external stratification only; "
            "not sample phenotype labels."
        ),
        "rows": matched.to_dict(orient="records") if not matched.empty else [],
    }
    json_path = report_dir / "study_atlas_enrichment.json"
    json_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    lines = [
        "# Study Atlas enrichment",
        "",
        f"- Catalog studies: **{payload['n_catalog_studies']}**",
        f"- Matched to Atlas: **{payload['n_matched']}**",
        "",
        "## Join methods",
        "",
    ]
    for method, count in sorted(by_method.items()):
        lines.append(f"- `{method}`: {count}")
    lines.extend(
        [
            "",
            "Atlas joins use curated GSE↔ES/PMID map and/or PMID bridge — never raw "
            "GSE = ES* equality.",
            "",
            "This is **external stratification metadata**, not training labels.",
            "",
        ]
    )
    if not matched.empty:
        lines.extend(
            [
                "| study_id | method | atlas_study_ids | pmid | cohorts | sample_size | tissues |",
                "| --- | --- | --- | --- | ---: | ---: | --- |",
            ]
        )
        preview = matched.head(50)
        for _, row in preview.iterrows():
            tissues = json.loads(str(row["tissues"]))
            tissue_preview = ", ".join(tissues[:2])
            if len(tissues) > 2:
                tissue_preview += ", …"
            lines.append(
                "| `{study_id}` | `{method}` | `{atlas}` | {pmid} | {cohorts} | {size} | {tissues} |".format(
                    study_id=row["study_id"],
                    method=row["join_method"],
                    atlas=str(row["atlas_study_ids"]),
                    pmid=row["pmid"] or "—",
                    cohorts=row["n_atlas_cohorts"] if row["n_atlas_cohorts"] is not None else "—",
                    size=row["total_sample_size"] if row["total_sample_size"] is not None else "—",
                    tissues=tissue_preview or "—",
                )
            )
        if len(matched) > len(preview):
            lines.append(f"| … | ({len(matched) - len(preview)} more) | | | | | |")
    md_path = report_dir / "study_atlas_enrichment.md"
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return json_path
