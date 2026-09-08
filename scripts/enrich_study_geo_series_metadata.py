#!/usr/bin/env python3
"""Enrich catalog studies with GEO series title/summary/design via NCBI + SOFT headers."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from mbs.geo_series_metadata import (
    build_geo_series_frame,
    fetch_gds_esummary,
    merge_geo_series_into_studies,
    overlay_soft_headers,
    write_geo_series_enrichment_report,
    write_geo_series_parquet,
)
from mbs.paths import DataPaths
from mbs.release import RELEASE_ID, release_inspection_report_dir, release_paths


def _catalog_gse_ids(study_parquet: Path) -> list[str]:
    if not study_parquet.is_file():
        raise FileNotFoundError(study_parquet)
    studies = pd.read_parquet(study_parquet, columns=["study_id"])
    return sorted(
        {
            str(s).strip().upper()
            for s in studies["study_id"].dropna().astype(str)
            if str(s).strip().upper().startswith("GSE")
        }
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional cap on GSE count (smoke / debug)",
    )
    parser.add_argument(
        "--skip-ncbi",
        action="store_true",
        help="Only parse cached SOFT headers (no NCBI esummary)",
    )
    parser.add_argument(
        "--merge-study-parquet",
        action="store_true",
        help="Also write merged metadata_json back into release study.parquet",
    )
    args = parser.parse_args()

    paths = DataPaths.from_environment()
    rp = release_paths(paths.data_root, RELEASE_ID)
    gse_ids = _catalog_gse_ids(rp.catalog_tables / "study.parquet")
    if args.limit is not None:
        gse_ids = gse_ids[: max(0, args.limit)]
    print(f"catalog GSE: {len(gse_ids)}", flush=True)

    rows: dict[str, dict] = {}
    if not args.skip_ncbi:
        rows = fetch_gds_esummary(gse_ids)
        print(f"ncbi esummary hits: {len(rows)}", flush=True)
    else:
        rows = {gse: {"study_id": gse, "source": "soft_header_only"} for gse in gse_ids}

    rows = overlay_soft_headers(rows, cache_root=paths.cache_root)
    n_design = sum(1 for r in rows.values() if r.get("overall_design"))
    print(f"with overall_design (SOFT): {n_design}", flush=True)

    frame = build_geo_series_frame(rows)
    out = write_geo_series_parquet(paths.data_root, frame)
    print(f"wrote {out}", flush=True)

    report_dir = release_inspection_report_dir(paths.project_root, RELEASE_ID)
    write_geo_series_enrichment_report(
        series_frame=frame,
        n_catalog_gse=len(gse_ids),
        report_dir=report_dir,
    )
    print(f"report {report_dir / 'geo_series_enrichment.md'}", flush=True)

    if args.merge_study_parquet and rp.catalog_tables.joinpath("study.parquet").is_file():
        studies = pd.read_parquet(rp.catalog_tables / "study.parquet")
        merged = merge_geo_series_into_studies(studies, frame)
        merged.to_parquet(rp.catalog_tables / "study.parquet", index=False)
        print(f"updated {rp.catalog_tables / 'study.parquet'}", flush=True)


if __name__ == "__main__":
    main()
