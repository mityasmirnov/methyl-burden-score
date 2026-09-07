#!/usr/bin/env python3
"""Fetch EWAS DataHub repository metadata census (paginated by platform)."""

from __future__ import annotations

import argparse

from mbs.datahub_census import (
    DATAHUB_PLATFORMS,
    fetch_datahub_repository_census,
)
from mbs.paths import DataPaths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--page-size", type=int, default=200)
    parser.add_argument("--delay-s", type=float, default=0.2)
    parser.add_argument("--from-cache-only", action="store_true")
    parser.add_argument(
        "--platforms",
        nargs="+",
        default=list(DATAHUB_PLATFORMS),
        help="Platform filters to page (default: 450K 850K 935K)",
    )
    args = parser.parse_args()

    paths = DataPaths.from_environment()
    cache_root = paths.cache_root / "ewas_datahub_repository"
    output_dir = paths.data_root / "canonical" / "phenotypes"
    manifest = fetch_datahub_repository_census(
        cache_root=cache_root,
        output_dir=output_dir,
        page_size=args.page_size,
        delay_s=args.delay_s,
        platforms=tuple(args.platforms),
        from_cache_only=args.from_cache_only,
    )
    # CLI status for operators (scripts are allowed to print).
    print(f"n_samples={manifest['n_samples']} official={manifest['official_n_samples']}")  # noqa: T201
    print(f"platform_totals_api={manifest['platform_totals_api']}")  # noqa: T201
    print(f"n_tissues={manifest['n_tissues']} n_diseases={manifest['n_diseases']}")  # noqa: T201
    print(f"wrote {manifest['census_parquet']}")  # noqa: T201


if __name__ == "__main__":
    main()
