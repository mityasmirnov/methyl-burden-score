"""GEO series-level metadata for study.metadata_json (not sample labels)."""

from __future__ import annotations

import gzip
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

import pandas as pd

from mbs.annotation.manifest import sha256_file, utc_now_iso, write_json
from mbs.atlas_study_enrichment import gds_uid_for_gse

NCBI_AGENT = "methyl-burden-score/1.0 (research pipeline)"
NCBI_DELAY_S = 0.34
SERIES_PARQUET_REL = Path("canonical/phenotypes/geo_series_metadata.parquet")
MAX_TEXT = 8000


def geo_series_parquet_path(data_root: Path) -> Path:
    return data_root / SERIES_PARQUET_REL


def _clip(text: object) -> str | None:
    if text is None or (isinstance(text, float) and pd.isna(text)):
        return None
    s = str(text).strip()
    if not s or s.lower() in {"nan", "none", "na", "<na>"}:
        return None
    if len(s) > MAX_TEXT:
        return s[: MAX_TEXT - 1] + "…"
    return s


def parse_series_soft_header(text: str) -> dict[str, Any]:
    """Parse Series_* fields from family SOFT; stop at first ^SAMPLE."""
    series: dict[str, list[str]] = {}
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith("^SAMPLE"):
            break
        if line.startswith("^SERIES") or line.startswith("^"):
            continue
        if not line.startswith("!"):
            continue
        if "=" not in line:
            continue
        key, _, value = line.partition("=")
        key = key.lstrip("!").strip()
        value = value.strip()
        if not key.startswith("Series_") and not key.startswith("series_"):
            # Accept bare keys after Series block starts
            if key.lower() in {
                "title",
                "summary",
                "overall_design",
                "pubmed_id",
                "geo_accession",
                "type",
                "platform_id",
            }:
                pass
            else:
                continue
        series.setdefault(key, []).append(value)

    def first(*keys: str) -> str | None:
        for key in keys:
            vals = series.get(key) or []
            for v in vals:
                clipped = _clip(v)
                if clipped:
                    return clipped
        return None

    pubmed_vals = series.get("Series_pubmed_id") or series.get("pubmed_id") or []
    pubmed_ids: list[str] = []
    for raw in pubmed_vals:
        pubmed_ids.extend(p.strip() for p in re.split(r"[;,]\s*", raw) if p.strip())
    platform_ids = [
        p.strip().upper()
        for p in (series.get("Series_platform_id") or series.get("platform_id") or [])
        if p.strip()
    ]
    study_id = first("Series_geo_accession", "geo_accession")
    return {
        "study_id": (study_id or "").strip().upper() or None,
        "title": first("Series_title", "title"),
        "summary": first("Series_summary", "summary"),
        "overall_design": first("Series_overall_design", "overall_design"),
        "series_type": first("Series_type", "type"),
        "pubmed_ids": sorted(dict.fromkeys(pubmed_ids)),
        "platform_ids": sorted(dict.fromkeys(platform_ids)),
        "source": "soft_header",
    }


def read_series_soft_header(path: Path) -> dict[str, Any]:
    """Read gzip family SOFT series header only (stream until first sample)."""
    lines: list[str] = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="latin-1", errors="replace") as handle:  # type: ignore[operator]
        for line in handle:
            lines.append(line)
            if line.startswith("^SAMPLE"):
                break
            # Cap pathological headers.
            if len(lines) > 50_000:
                break
    return parse_series_soft_header("".join(lines))


def _ncbi_get(url: str) -> dict[str, Any]:
    req = urllib.request.Request(url, headers={"User-Agent": NCBI_AGENT})  # noqa: S310
    with urllib.request.urlopen(req, timeout=60) as resp:  # noqa: S310
        return json.loads(resp.read().decode("utf-8"))


def fetch_gds_esummary(
    gse_ids: list[str],
    *,
    batch_size: int = 100,
    delay_s: float = NCBI_DELAY_S,
) -> dict[str, dict[str, Any]]:
    """Fetch NCBI GDS esummary for GSE accessions → series context fields."""
    out: dict[str, dict[str, Any]] = {}
    ordered = sorted({g.strip().upper() for g in gse_ids if g.strip().upper().startswith("GSE")})
    for start in range(0, len(ordered), batch_size):
        batch = ordered[start : start + batch_size]
        uids = [gds_uid_for_gse(gse) for gse in batch]
        params = urllib.parse.urlencode({"db": "gds", "retmode": "json", "id": ",".join(uids)})
        payload = _ncbi_get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?{params}")
        result = payload.get("result") or {}
        for gse, uid in zip(batch, uids, strict=True):
            record = result.get(uid) or {}
            accession = str(record.get("accession") or "").strip().upper()
            if accession != gse:
                continue
            gpl_raw = record.get("gpl")
            platform_ids: list[str] = []
            if isinstance(gpl_raw, list):
                platform_ids = [f"GPL{str(x).strip()}" for x in gpl_raw if str(x).strip()]
            elif gpl_raw is not None and str(gpl_raw).strip():
                # Sometimes a single id or semicolon list.
                for part in re.split(r"[;,]\s*", str(gpl_raw)):
                    part = part.strip()
                    if not part:
                        continue
                    platform_ids.append(part if part.upper().startswith("GPL") else f"GPL{part}")
            pmids = [str(p).strip() for p in (record.get("pubmedids") or []) if str(p).strip()]
            out[gse] = {
                "study_id": gse,
                "title": _clip(record.get("title") or record.get("seriestitle")),
                "summary": _clip(record.get("summary")),
                "overall_design": None,
                "series_type": _clip(record.get("gdstype")),
                "pubmed_ids": sorted(dict.fromkeys(pmids)),
                "platform_ids": sorted(dict.fromkeys(p.upper() for p in platform_ids)),
                "n_samples_ncbi": record.get("n_samples"),
                "taxon": _clip(record.get("taxon") or record.get("samplestaxa")),
                "platform_title": _clip(record.get("platformtitle")),
                "source": "ncbi_esummary",
            }
        if start + batch_size < len(ordered):
            time.sleep(delay_s)
    return out


def overlay_soft_headers(
    rows: dict[str, dict[str, Any]],
    *,
    cache_root: Path,
) -> dict[str, dict[str, Any]]:
    """Fill overall_design / missing fields from cached family SOFT headers."""
    out = {k: dict(v) for k, v in rows.items()}
    for gse, row in list(out.items()):
        soft = cache_root / "geo" / gse / f"{gse}_family.soft.gz"
        if not soft.is_file():
            continue
        try:
            header = read_series_soft_header(soft)
        except OSError:
            continue
        if header.get("overall_design"):
            row["overall_design"] = header["overall_design"]
        if not row.get("title") and header.get("title"):
            row["title"] = header["title"]
        if not row.get("summary") and header.get("summary"):
            row["summary"] = header["summary"]
        if not row.get("series_type") and header.get("series_type"):
            row["series_type"] = header["series_type"]
        if header.get("pubmed_ids"):
            row["pubmed_ids"] = sorted(
                dict.fromkeys([*(row.get("pubmed_ids") or []), *header["pubmed_ids"]])
            )
        if header.get("platform_ids"):
            row["platform_ids"] = sorted(
                dict.fromkeys([*(row.get("platform_ids") or []), *header["platform_ids"]])
            )
        sources = {row.get("source"), "soft_header"}
        row["source"] = "+".join(sorted(s for s in sources if s))
        out[gse] = row
    return out


def build_geo_series_frame(rows: dict[str, dict[str, Any]]) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    for gse, row in sorted(rows.items()):
        records.append(
            {
                "study_id": gse,
                "title": row.get("title"),
                "summary": row.get("summary"),
                "overall_design": row.get("overall_design"),
                "series_type": row.get("series_type"),
                "pubmed_ids": json.dumps(row.get("pubmed_ids") or [], sort_keys=True),
                "platform_ids": json.dumps(row.get("platform_ids") or [], sort_keys=True),
                "n_samples_ncbi": row.get("n_samples_ncbi"),
                "taxon": row.get("taxon"),
                "platform_title": row.get("platform_title"),
                "source": row.get("source"),
                "fetched_at": utc_now_iso(),
            }
        )
    return pd.DataFrame(records)


def write_geo_series_parquet(data_root: Path, frame: pd.DataFrame) -> Path:
    out = geo_series_parquet_path(data_root)
    out.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(out, index=False)
    write_json(
        out.with_suffix(".manifest.json"),
        {
            "path": str(out),
            "sha256": sha256_file(out),
            "n_studies": int(len(frame)),
            "created_at": utc_now_iso(),
        },
    )
    return out


def load_geo_series_frame(data_root: Path) -> pd.DataFrame:
    path = geo_series_parquet_path(data_root)
    if not path.is_file():
        return pd.DataFrame()
    return pd.read_parquet(path)


def merge_geo_series_into_studies(
    studies: pd.DataFrame,
    series_frame: pd.DataFrame,
) -> pd.DataFrame:
    """Merge series context under study.metadata_json.geo (non-destructive)."""
    if studies.empty or series_frame.empty:
        return studies
    out = studies.copy()
    if "metadata_json" not in out.columns:
        out["metadata_json"] = None
    by_id = series_frame.set_index("study_id", drop=False)
    for idx, row in out.iterrows():
        study_id = str(row.get("study_id") or "").strip().upper()
        if study_id not in by_id.index:
            continue
        ser = by_id.loc[study_id]
        if isinstance(ser, pd.DataFrame):
            ser = ser.iloc[0]
        meta_raw = row.get("metadata_json")
        meta: dict[str, Any] = {}
        if meta_raw is not None and not (isinstance(meta_raw, float) and pd.isna(meta_raw)):
            try:
                parsed = json.loads(str(meta_raw))
                if isinstance(parsed, dict):
                    meta = parsed
            except json.JSONDecodeError:
                meta = {}
        geo = dict(meta.get("geo") or {})
        for key in (
            "title",
            "summary",
            "overall_design",
            "series_type",
            "taxon",
            "platform_title",
            "n_samples_ncbi",
            "source",
        ):
            val = ser.get(key)
            if val is None or (isinstance(val, float) and pd.isna(val)):
                continue
            if key == "n_samples_ncbi":
                try:
                    geo[key] = int(val)
                except (TypeError, ValueError):
                    continue
            else:
                text = _clip(val)
                if text is not None and text != "":
                    geo[key] = text
        for list_key in ("pubmed_ids", "platform_ids"):
            raw = ser.get(list_key)
            values: list[str] = []
            if isinstance(raw, str) and raw.strip():
                try:
                    loaded = json.loads(raw)
                    if isinstance(loaded, list):
                        values = [str(v) for v in loaded if str(v).strip()]
                except json.JSONDecodeError:
                    values = []
            if values:
                existing = geo.get(list_key) or []
                if not isinstance(existing, list):
                    existing = []
                geo[list_key] = sorted(dict.fromkeys([*existing, *values]))
        fetched = ser.get("fetched_at")
        if fetched is not None and not (isinstance(fetched, float) and pd.isna(fetched)):
            geo["series_fetched_at"] = str(fetched)
        meta["geo"] = geo
        out.at[idx, "metadata_json"] = json.dumps(meta, sort_keys=True)
    return out


def write_geo_series_enrichment_report(
    *,
    series_frame: pd.DataFrame,
    n_catalog_gse: int,
    report_dir: Path,
) -> Path:
    report_dir.mkdir(parents=True, exist_ok=True)
    n = int(len(series_frame))
    n_title = int(series_frame["title"].notna().sum()) if n else 0
    n_summary = int(series_frame["summary"].notna().sum()) if n else 0
    n_design = int(series_frame["overall_design"].notna().sum()) if n else 0
    payload = {
        "generated_at": utc_now_iso(),
        "n_catalog_gse": n_catalog_gse,
        "n_series_enriched": n,
        "n_with_title": n_title,
        "n_with_summary": n_summary,
        "n_with_overall_design": n_design,
        "sources": (
            series_frame["source"].value_counts(dropna=False).to_dict() if n else {}
        ),
    }
    write_json(report_dir / "geo_series_enrichment.json", payload)
    lines = [
        "# GEO series enrichment (study context)",
        "",
        f"- Generated: `{payload['generated_at']}`",
        f"- Catalog GSE studies targeted: **{n_catalog_gse}**",
        f"- Series rows written: **{n}**",
        f"- With title: **{n_title}**",
        f"- With summary: **{n_summary}**",
        f"- With overall_design (SOFT header): **{n_design}**",
        "",
        "Stored under `study.metadata_json.geo` (title / summary / overall_design / "
        "series_type / pubmed_ids). Not sample phenotype labels; not encoder features.",
        "",
        "## Sources",
        "",
    ]
    for src, count in sorted((payload["sources"] or {}).items(), key=lambda x: (-int(x[1]), str(x[0]))):
        lines.append(f"- `{src}`: {count}")
    lines.append("")
    (report_dir / "geo_series_enrichment.md").write_text("\n".join(lines), encoding="utf-8")
    return report_dir
