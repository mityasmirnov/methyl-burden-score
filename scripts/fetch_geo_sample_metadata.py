#!/usr/bin/env python3
"""Fetch GEO family SOFT for pilot GSE list → geo_sample_metadata.parquet."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from mbs.annotation.manifest import sha256_file, utc_now_iso, write_json
from mbs.geo_metadata import (
    build_geo_frame_from_soft,
    consolidate_geo_sample_rows,
    download_family_soft,
    load_geo_tissue_aliases,
    read_cached_soft,
    resolve_tissue_ontology_path,
    write_geo_parquet,
)
from mbs.paths import DataPaths
from mbs.training.phenotype_table import load_tissue_ontology

DEFAULT_STUDIES = Path("configs/data/geo_backfill_pilot_gse.txt")


def _report_subdir(studies_file: Path) -> str:
    name = studies_file.name.lower()
    if "batch" in name:
        return "geo_backfill_batch"
    return "geo_backfill_pilot"


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
    parser.add_argument(
        "--tissue-ontology",
        type=Path,
        default=None,
        help="Tissue ontology YAML (default: Hub nine-pack under data root)",
    )
    parser.add_argument(
        "--tissue-aliases",
        type=Path,
        default=None,
        help="GEO tissue alias YAML (default: configs/data/geo_tissue_aliases.yaml)",
    )
    args = parser.parse_args()

    paths = DataPaths.from_environment()
    paths.ensure_directories()
    study_ids = _load_study_ids(args.studies_file.resolve())
    if not study_ids:
        raise SystemExit(f"no studies in {args.studies_file}")

    ont_path = args.tissue_ontology or resolve_tissue_ontology_path(paths.data_root)
    ontology = None
    if ont_path is not None and ont_path.is_file():
        ontology = load_tissue_ontology(ont_path)
    aliases = load_geo_tissue_aliases(args.tissue_aliases)
    if ontology is None:
        sys.stdout.write("warn: no tissue ontology found; tissues stay unmapped\n")

    frames: list[pd.DataFrame] = []
    failures: list[dict[str, str]] = []
    per_study: list[dict[str, object]] = []
    fetched_at = utc_now_iso()

    for gse in study_ids:
        cache_path = paths.cache_root / "geo" / gse / f"{gse}_family.soft.gz"
        try:
            if args.from_cache_only:
                if not cache_path.is_file():
                    failures.append({"study_id": gse, "error": "cache missing"})
                    per_study.append(
                        {"study_id": gse, "download_status": "fail", "error": "cache missing"}
                    )
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
            frame = build_geo_frame_from_soft(
                text,
                fetched_at=fetched_at,
                soft_sha256=digest,
                ontology=ontology,
                aliases=aliases,
            )
            tissue_stats = dict(frame.attrs.get("tissue_map_stats") or {})
            pheno_counts = {
                "age": int(frame["age"].notna().sum()) if "age" in frame.columns else 0,
                "sex": int(frame["sex"].notna().sum()) if "sex" in frame.columns else 0,
                "tissue": int(
                    (frame["tissue_map_status"] == "mapped").sum()
                    if "tissue_map_status" in frame.columns
                    else 0
                ),
                "disease": int(frame["disease"].notna().sum())
                if "disease" in frame.columns
                else 0,
                "cancer": int(frame["cancer"].notna().sum()) if "cancer" in frame.columns else 0,
            }
            per_study.append(
                {
                    "study_id": gse,
                    "download_status": "ok",
                    "error": None,
                    "soft_sha256": digest,
                    "n_geo_gsm": len(frame),
                    "tissue_map": tissue_stats,
                    "phenotype_counts": pheno_counts,
                }
            )
            if not frame.empty:
                frames.append(frame)
            sys.stdout.write(f"ok {gse} n_samples={len(frame)}\n")
        except OSError as exc:
            failures.append({"study_id": gse, "error": str(exc)})
            per_study.append({"study_id": gse, "download_status": "fail", "error": str(exc)})
            sys.stdout.write(f"fail {gse}: {exc}\n")

    if not frames:
        raise SystemExit("no GEO samples fetched")

    combined = pd.concat(frames, ignore_index=True)
    combined, conflict_stats = consolidate_geo_sample_rows(combined)
    out = write_geo_parquet(paths.data_root, combined, study_ids=study_ids)
    status_path = (
        paths.project_root
        / "reports"
        / "inspection"
        / "deepmat_data_v1"
        / _report_subdir(args.studies_file)
        / "fetch_status.json"
    )
    status_path.parent.mkdir(parents=True, exist_ok=True)
    write_json(
        status_path,
        {
            "generated_at": fetched_at,
            "n_samples": len(combined),
            "tissue_ontology": str(ont_path) if ont_path else None,
            "conflict_stats": {
                k: v for k, v in conflict_stats.items() if k != "conflicts"
            },
            "conflicts": conflict_stats.get("conflicts") or [],
            "per_study": per_study,
            "failures": failures,
        },
    )
    sys.stdout.write(
        f"wrote {out} n_samples={len(combined)} "
        f"conflicts={conflict_stats.get('n_conflict_samples', 0)} "
        f"failures={len(failures)}\n"
    )
    if failures:
        for item in failures:
            sys.stdout.write(f"  {item['study_id']}: {item['error']}\n")
    if conflict_stats.get("n_conflict_samples"):
        sys.stdout.write(
            f"  metadata conflicts: {conflict_stats['n_conflict_samples']} GSM "
            f"({conflict_stats['n_conflict_fields']} fields); see {status_path}\n"
        )


if __name__ == "__main__":
    main()
