#!/usr/bin/env python3
"""Fetch GEO family SOFT for pilot GSE list → geo_sample_metadata.parquet."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from mbs.annotation.manifest import sha256_file, utc_now_iso
from mbs.geo_metadata import (
    build_geo_frame_from_soft,
    download_family_soft,
    read_cached_soft,
    write_geo_parquet,
)
from mbs.paths import DataPaths

DEFAULT_STUDIES = Path("configs/data/geo_backfill_pilot_gse.txt")


def _load_study_ids(path: Path) -> list[str]:
    lines = path.read_text(encoding="utf-8").splitlines()
    return [line.strip().upper() for line in lines if line.strip() and not line.startswith("#")]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--studies-file",
        type=Path,
        default=DEFAULT_STUDIES,
        help="Newline-delimited GSE accessions (default: pilot list)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-download cached SOFT even if present",
    )
    parser.add_argument(
        "--from-cache-only",
        action="store_true",
        help="Use cached SOFT only; do not download missing studies",
    )
    args = parser.parse_args()

    paths = DataPaths.from_environment()
    paths.ensure_directories()
    study_ids = _load_study_ids(args.studies_file.resolve())
    if not study_ids:
        raise SystemExit(f"no studies in {args.studies_file}")

    frames: list[pd.DataFrame] = []
    failures: list[dict[str, str]] = []
    fetched_at = utc_now_iso()

    for gse in study_ids:
        cache_path = paths.cache_root / "geo" / gse / f"{gse}_family.soft.gz"
        try:
            if args.from_cache_only:
                if not cache_path.is_file():
                    failures.append({"study_id": gse, "error": "cache missing"})
                    continue
                digest = sha256_file(cache_path)
                text = read_cached_soft(cache_path)
            else:
                cache_path, digest = download_family_soft(
                    gse,
                    cache_root=paths.cache_root,
                    force=args.force,
                )
                text = read_cached_soft(cache_path)
            frame = build_geo_frame_from_soft(text, fetched_at=fetched_at, soft_sha256=digest)
            if not frame.empty:
                frames.append(frame)
            sys.stdout.write(f"ok {gse} n_samples={len(frame)}\n")
        except OSError as exc:
            failures.append({"study_id": gse, "error": str(exc)})
            sys.stdout.write(f"fail {gse}: {exc}\n")

    if not frames:
        raise SystemExit("no GEO samples fetched")

    combined = pd.concat(frames, ignore_index=True)
    combined = combined.drop_duplicates(subset=["sample_id"], keep="first")
    out = write_geo_parquet(paths.data_root, combined, study_ids=study_ids)
    sys.stdout.write(f"wrote {out} n_samples={len(combined)} failures={len(failures)}\n")
    if failures:
        for item in failures:
            sys.stdout.write(f"  {item['study_id']}: {item['error']}\n")


if __name__ == "__main__":
    main()
