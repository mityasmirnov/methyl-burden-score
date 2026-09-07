#!/usr/bin/env python3
"""Remap GEO tissue labels from source_name + aliases without re-parsing SOFT.

Use after expanding ``configs/data/geo_tissue_aliases.yaml`` or when parquet was
built before the source_name fallback. Avoids loading all family SOFT into RAM.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from mbs.annotation.manifest import sha256_file, utc_now_iso, write_json
from mbs.geo_metadata import (
    load_geo_tissue_aliases,
    remap_geo_tissue_frame,
    resolve_tissue_ontology_path,
)
from mbs.paths import DataPaths
from mbs.training.phenotype_table import load_tissue_ontology


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--parquet",
        type=Path,
        default=None,
        help="geo_sample_metadata.parquet (default under data/canonical/phenotypes/)",
    )
    parser.add_argument(
        "--tissue-aliases",
        type=Path,
        default=None,
        help="Alias YAML (default: configs/data/geo_tissue_aliases.yaml)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Remap every row (apply alias/ontology upgrades)",
    )
    args = parser.parse_args()
    paths = DataPaths.from_environment()
    parquet = args.parquet or (paths.data_root / "canonical/phenotypes/geo_sample_metadata.parquet")
    if not parquet.is_file():
        raise SystemExit(f"missing {parquet}")
    ont_path = resolve_tissue_ontology_path(paths.data_root)
    if ont_path is None or not ont_path.is_file():
        raise SystemExit("tissue ontology not found")
    ontology = load_tissue_ontology(ont_path)
    aliases = load_geo_tissue_aliases(
        args.tissue_aliases or paths.project_root / "configs/data/geo_tissue_aliases.yaml"
    )
    df = pd.read_parquet(parquet)
    out, stats = remap_geo_tissue_frame(
        df, ontology=ontology, aliases=aliases, force=args.force
    )
    out.to_parquet(parquet, index=False)
    write_json(
        parquet.with_suffix(".manifest.json"),
        {
            "path": str(parquet),
            "sha256": sha256_file(parquet),
            "remapped_at": utc_now_iso(),
            **stats,
        },
    )
    print(f"before {stats['tissue_map_status_before']}", flush=True)
    print(f"after  {stats['tissue_map_status_after']}", flush=True)
    print(f"from_source_name {stats['n_from_source_name']}", flush=True)
    print(f"wrote {parquet}", flush=True)


if __name__ == "__main__":
    main()
