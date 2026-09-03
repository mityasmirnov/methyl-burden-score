#!/usr/bin/env python3
"""Seed curated GSE↔Atlas map from NCBI GEO PubMed IDs (verified PMID bridge)."""

from __future__ import annotations

import argparse
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

from mbs.atlas_study_enrichment import (
    build_atlas_reference,
    gds_uid_for_gse,
    load_atlas_gse_map,
)
from mbs.paths import DataPaths

NCBI_AGENT = "methyl-burden-score/1.0 (research pipeline)"
DEFAULT_MAP = Path("configs/data/atlas_gse_es_map.tsv")
DEFAULT_ATLAS = Path("data/raw/ewas_atlas")
DEFAULT_STUDY_PARQUET = Path(
    "data/canonical/releases/deepmat-data-v1/catalog/tables/study.parquet"
)


def _ncbi_get(url: str) -> dict:
    req = urllib.request.Request(url, headers={"User-Agent": NCBI_AGENT})
    with urllib.request.urlopen(req, timeout=60) as resp:
        return json.loads(resp.read().decode("utf-8"))


def fetch_geo_pubmed_ids(gse_ids: list[str], *, batch_size: int = 200, delay_s: float = 0.34) -> dict[str, list[str]]:
    """Return GSE → PubMed ID list from NCBI GEO DataSets esummary."""
    out: dict[str, list[str]] = {}
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
            pmids = [str(p).strip() for p in record.get("pubmedids") or [] if str(p).strip()]
            if pmids:
                out[gse] = sorted(dict.fromkeys(pmids))
        if start + batch_size < len(ordered):
            time.sleep(delay_s)
    return out


def build_map_from_geo(
    *,
    gse_ids: list[str],
    atlas_root: Path,
    geo_pmids: dict[str, list[str]],
) -> pd.DataFrame:
    """Build map rows where GEO PubMed IDs hit Atlas studies."""
    _, _, pmid_to_study = build_atlas_reference(atlas_root=atlas_root)
    rows: list[dict[str, str]] = []
    for gse in sorted(set(gse_ids)):
        gse = gse.strip().upper()
        if not gse.startswith("GSE"):
            continue
        pmids = geo_pmids.get(gse) or []
        atlas_ids: list[str] = []
        matched_pmids: list[str] = []
        for pmid in pmids:
            hits = pmid_to_study.get(pmid) or []
            if hits:
                matched_pmids.append(pmid)
                atlas_ids.extend(hits)
        if not matched_pmids:
            continue
        rows.append(
            {
                "gse_id": gse,
                "atlas_study_id": ";".join(sorted(dict.fromkeys(atlas_ids))),
                "pmid": matched_pmids[0] if len(matched_pmids) == 1 else ";".join(matched_pmids),
                "source": "geo_ncbi",
            }
        )
    return pd.DataFrame(rows)


def merge_map(existing: pd.DataFrame, seeded: pd.DataFrame) -> pd.DataFrame:
    """Keep manual rows; replace geo_ncbi rows; add new seeds."""
    columns = ["gse_id", "atlas_study_id", "pmid", "source"]
    if existing.empty:
        merged = seeded.copy()
    else:
        manual = existing.loc[existing["source"] != "geo_ncbi"].copy()
        manual_gse = set(manual["gse_id"].astype(str))
        seeded = seeded.loc[~seeded["gse_id"].isin(manual_gse)].copy()
        merged = pd.concat([manual, seeded], ignore_index=True)
    merged = merged[columns].sort_values(by=["gse_id"]).reset_index(drop=True)
    return merged


def write_map(path: Path, frame: pd.DataFrame) -> None:
    header = (
        "# Optional curated GSE↔Atlas bridge for study-level enrichment.\n"
        "# Never join Hub project_id (GSE*) to Atlas study_ID (ES*) by raw string equality.\n"
        "# Rows with source=geo_ncbi are seeded from NCBI GEO pubmedids + Atlas PMID index.\n"
        "#\n"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    body = frame.to_csv(sep="\t", index=False)
    path.write_text(header + body, encoding="utf-8")


def catalog_gse_ids(study_parquet: Path) -> list[str]:
    if not study_parquet.is_file():
        raise FileNotFoundError(study_parquet)
    studies = pd.read_parquet(study_parquet, columns=["study_id"])
    return [
        str(s)
        for s in studies["study_id"].dropna().astype(str)
        if str(s).upper().startswith("GSE")
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--map-path", type=Path, default=DEFAULT_MAP)
    parser.add_argument("--atlas-root", type=Path, default=DEFAULT_ATLAS)
    parser.add_argument("--study-parquet", type=Path, default=DEFAULT_STUDY_PARQUET)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    paths = DataPaths.from_environment()
    map_path = args.map_path if args.map_path.is_absolute() else paths.project_root / args.map_path
    atlas_root = args.atlas_root if args.atlas_root.is_absolute() else paths.data_root / "raw" / "ewas_atlas"
    study_parquet = (
        args.study_parquet
        if args.study_parquet.is_absolute()
        else paths.data_root / "canonical/releases/deepmat-data-v1/catalog/tables/study.parquet"
    )

    gse_ids = catalog_gse_ids(study_parquet)
    existing = load_atlas_gse_map(map_path)
    geo_pmids = fetch_geo_pubmed_ids(gse_ids)
    seeded = build_map_from_geo(gse_ids=gse_ids, atlas_root=atlas_root, geo_pmids=geo_pmids)
    merged = merge_map(existing, seeded)

    print(f"catalog GSE studies: {len(gse_ids)}")
    print(f"GEO with PubMed IDs: {len(geo_pmids)}")
    print(f"seeded Atlas matches: {len(seeded)}")
    print(f"map rows (total): {len(merged)}")

    if args.dry_run:
        print(merged.head(10).to_string(index=False))
        return

    write_map(map_path, merged)
    print(f"wrote {map_path}")


if __name__ == "__main__":
    main()
