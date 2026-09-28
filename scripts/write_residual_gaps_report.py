#!/usr/bin/env python3
"""Write live residual catalog / DataHub / platform gap report.

Computes census∩catalog gaps, empty EWAS_db dirs (with optional remote probe
cache), and null study.platform_id classification. Does not download.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import duckdb
import pandas as pd

from mbs.annotation.manifest import utc_now_iso, write_json
from mbs.datahub_census import load_datahub_census
from mbs.paths import DataPaths
from mbs.platform_id import normalize_platform
from mbs.sample_overview_enrich import (
    fill_null_study_platforms_from_ewas_db,
    infer_study_platform_from_assay_dir,
)

RELEASE_ID = "deepmat-data-v1"


def _blank(value: object) -> bool:
    if value is None:
        return True
    if isinstance(value, float) and pd.isna(value):
        return True
    text = str(value).strip()
    return text == "" or text.lower() == "nan"


def _explode_census_sample_ids(census: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    for rec in census.to_dict(orient="records"):
        study = str(rec.get("project_id", "")).strip()
        raw = str(rec.get("sample_id", "")).strip()
        plat = rec.get("catalog_platform_id")
        for part in (p.strip() for p in raw.split(",") if p.strip()):
            rows.append(
                {
                    "sample_id": part,
                    "study_id": study,
                    "catalog_platform_id": (
                        str(plat) if plat is not None and not _blank(plat) else ""
                    ),
                }
            )
    return pd.DataFrame(rows)


def _classify_null_platforms(
    *,
    null_studies: pd.DataFrame,
    census: pd.DataFrame,
    ewas_db_root: Path,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for study_id, n_samples in null_studies.itertuples(index=False):
        sub = census[census["project_id"].astype(str) == str(study_id)]
        plats = sorted(
            {
                str(p).strip()
                for p in sub.get("catalog_platform_id", pd.Series(dtype=object)).dropna()
                if not _blank(p)
            }
        )
        fill: str | None = None
        if len(plats) > 1:
            reason = "mixed_census"
        elif len(plats) == 1:
            reason = "single_census_unapplied"
            fill = plats[0]
        else:
            inferred = infer_study_platform_from_assay_dir(ewas_db_root / str(study_id))
            if inferred is not None:
                reason = "assay_probe_count"
                fill = inferred
            else:
                reason = "no_single_platform_source"
        rows.append(
            {
                "study_id": str(study_id),
                "n_samples": int(n_samples),
                "reason": reason,
                "fill": fill,
                "census_platforms": "|".join(plats),
            }
        )
    return pd.DataFrame(rows)


def build_residual_gaps(
    *,
    paths: DataPaths,
    release_id: str = RELEASE_ID,
    apply_platform_fill: bool = False,
) -> dict[str, Any]:
    rp_catalog = (
        paths.data_root / "canonical" / "releases" / release_id / "catalog" / "catalog.duckdb"
    )
    ewas_root = paths.data_root / "raw" / "ewas_datahub" / "EWAS_db"
    study_parquet = (
        paths.data_root
        / "canonical"
        / "releases"
        / release_id
        / "catalog"
        / "tables"
        / "study.parquet"
    )

    con = duckdb.connect(str(rp_catalog), read_only=not apply_platform_fill)
    n_samples = int(con.execute("select count(*) from sample").fetchone()[0])
    n_studies = int(con.execute("select count(*) from study").fetchone()[0])
    n_assay = int(con.execute("select count(*) from assay_file").fetchone()[0])
    plat_counts = con.execute(
        """
        select coalesce(platform_id, '(null)') as platform_id, count(*) as n
        from study group by 1 order by 2 desc
        """
    ).df()
    null_studies = con.execute(
        """
        select s.study_id, count(sa.sample_id) as n_samples
        from study s
        left join sample sa using (study_id)
        where s.platform_id is null
           or cast(s.platform_id as varchar) in ('', 'nan', 'NaN')
        group by 1
        order by 2 desc
        """
    ).df()

    census = load_datahub_census(paths.data_root)
    exploded = _explode_census_sample_ids(census)
    census_ids = set(exploded["sample_id"].astype(str))
    cat_ids = set(
        con.execute("select sample_id::varchar as sample_id from sample")
        .df()["sample_id"]
        .astype(str)
    )
    census_not_catalog = sorted(census_ids - cat_ids)
    catalog_not_census = sorted(cat_ids - census_ids)

    empty_dirs = sorted(
        p.name for p in ewas_root.iterdir() if p.is_dir() and not any(p.glob("*.txt"))
    )
    with_txt = sum(1 for p in ewas_root.iterdir() if p.is_dir() and any(p.glob("*.txt")))
    n_dirs = sum(1 for p in ewas_root.iterdir() if p.is_dir())
    n_txt = sum(1 for p in ewas_root.rglob("*.txt") if p.is_file())
    n_gsm = sum(1 for p in ewas_root.rglob("GSM*.txt") if p.is_file())

    probe_path = (
        paths.project_root
        / "reports"
        / "inspection"
        / "deepmat_data_v1"
        / "ewas_db_empty_remote_probe.json"
    )
    remote_empty: list[str] = []
    recoverable: list[dict[str, Any]] = []
    if probe_path.is_file():
        probe = json.loads(probe_path.read_text(encoding="utf-8"))
        remote_empty = list(probe.get("empty_both") or [])
        recoverable = list(probe.get("recoverable") or [])

    classification = _classify_null_platforms(
        null_studies=null_studies,
        census=census,
        ewas_db_root=ewas_root,
    )
    fillable = classification[
        classification["fill"].notna()
        & classification["reason"].isin(["assay_probe_count", "single_census_unapplied"])
    ]

    platform_fill_stats: dict[str, int] = {
        "n_studies_considered": 0,
        "n_studies_platform_set": 0,
    }
    if apply_platform_fill and not fillable.empty:
        allow = set(fillable["study_id"].astype(str))
        if study_parquet.is_file():
            studies = pd.read_parquet(study_parquet)
        else:
            studies = con.execute("select * from study").df()
        studies, platform_fill_stats = fill_null_study_platforms_from_ewas_db(
            studies,
            ewas_db_root=ewas_root,
            study_ids=allow,
        )
        study_parquet.parent.mkdir(parents=True, exist_ok=True)
        studies.to_parquet(study_parquet, index=False)
        # DuckDB parent-table UPDATE hits FK rewrite limits; parquet is SoT until
        # the next `mbs catalog refresh-release` (which reloads tables).
        platform_fill_stats["duckdb_update"] = "skipped_fk_limitation_use_refresh_release"
        sample_counts = {
            str(r.study_id): int(r.n_samples)
            for r in null_studies.itertuples(index=False)
        }
        null_mask = studies["platform_id"].map(_blank)
        null_studies = pd.DataFrame(
            [
                {
                    "study_id": sid,
                    "n_samples": sample_counts.get(sid, 0),
                }
                for sid in studies.loc[null_mask, "study_id"].astype(str)
            ]
        )
        plat_counts = (
            studies.assign(
                platform_id=studies["platform_id"].map(
                    lambda v: "(null)" if _blank(v) else str(v)
                )
            )
            .groupby("platform_id", as_index=False)
            .size()
            .rename(columns={"size": "n"})
            .sort_values("n", ascending=False)
        )
        classification = _classify_null_platforms(
            null_studies=null_studies,
            census=census,
            ewas_db_root=ewas_root,
        )

    con.close()

    miss_by_study = (
        exploded[exploded["sample_id"].isin(census_not_catalog)]
        .groupby("study_id")
        .size()
        .sort_values(ascending=False)
    )

    summary: dict[str, Any] = {
        "generated_at": utc_now_iso(),
        "release_id": release_id,
        "catalog": {
            "n_samples": n_samples,
            "n_studies": n_studies,
            "n_assay_files": n_assay,
            "platform_study_counts": {
                str(r.platform_id): int(r.n) for r in plat_counts.itertuples(index=False)
            },
        },
        "disk_ewas_db": {
            "n_dirs": n_dirs,
            "n_with_txt": with_txt,
            "n_empty": len(empty_dirs),
            "n_txt_files": n_txt,
            "n_gsm_files": n_gsm,
            "empty_dirs": empty_dirs,
        },
        "vs_datahub_census": {
            "n_census_unique_samples": len(census_ids),
            "n_catalog_samples": len(cat_ids),
            "overlap": len(census_ids & cat_ids),
            "census_not_in_catalog": len(census_not_catalog),
            "catalog_not_in_census": len(catalog_not_census),
            "census_not_in_catalog_top_studies": {
                str(k): int(v) for k, v in miss_by_study.head(20).items()
            },
        },
        "remote_probe": {
            "path": str(probe_path) if probe_path.is_file() else None,
            "n_empty_both_sides": len(remote_empty),
            "n_recoverable_delta_gt0": len(recoverable),
            "note": (
                "Empty local dirs with remote_txt=0 cannot be filled from the Hub mirror. "
                "Partial census>disk studies are Hub-side incomplete, not local download debt."
            ),
        },
        "platform_null": {
            "n_studies": int(len(null_studies)),
            "n_samples": int(null_studies["n_samples"].sum()) if not null_studies.empty else 0,
            "reason_counts": {
                str(k): int(v) for k, v in classification["reason"].value_counts().items()
            },
            "classification": classification.to_dict(orient="records"),
        },
        "assay_platform_fill": platform_fill_stats,
        "applied_platform_fill": apply_platform_fill,
    }
    return summary


def write_markdown(summary: dict[str, Any], path: Path) -> None:
    cat = summary["catalog"]
    disk = summary["disk_ewas_db"]
    vs = summary["vs_datahub_census"]
    plat = summary["platform_null"]
    remote = summary["remote_probe"]
    fill = summary.get("assay_platform_fill") or {}
    lines = [
        "# Residual catalog / DataHub / platform gaps",
        "",
        f"- Generated: `{summary['generated_at']}`",
        f"- Release: `{summary['release_id']}`",
        f"- Applied assay platform fill this run: **{summary.get('applied_platform_fill')}** "
        f"(set {fill.get('n_studies_platform_set', 0)} / considered {fill.get('n_studies_considered', 0)})",
        "",
        "## Catalog (release DuckDB)",
        "",
        f"- Samples: **{cat['n_samples']:,}**",
        f"- Studies: **{cat['n_studies']:,}**",
        f"- Assay files: **{cat['n_assay_files']:,}**",
        "",
        "| Platform | Studies |",
        "| --- | ---: |",
    ]
    for pid, n in cat["platform_study_counts"].items():
        lines.append(f"| `{pid}` | {n} |")
    lines.extend(
        [
            "",
            "## EWAS_db on disk",
            "",
            f"- Study dirs: **{disk['n_dirs']}** (with `*.txt`: **{disk['n_with_txt']}**, "
            f"empty: **{disk['n_empty']}**)",
            f"- Sample files: **{disk['n_txt_files']:,}** (`GSM*.txt`: **{disk['n_gsm_files']:,}**)",
            "",
            "## Catalog vs DataHub census",
            "",
            f"- Census unique sample IDs (comma-split): **{vs['n_census_unique_samples']:,}**",
            f"- Catalog samples: **{vs['n_catalog_samples']:,}**",
            f"- Overlap: **{vs['overlap']:,}**",
            f"- In census, not in catalog: **{vs['census_not_in_catalog']:,}**",
            f"- In catalog, not in census: **{vs['catalog_not_in_census']:,}**",
            "",
            "Top studies for census-not-in-catalog:",
            "",
            "| Study | Missing census IDs |",
            "| --- | ---: |",
        ]
    )
    for study, n in vs["census_not_in_catalog_top_studies"].items():
        lines.append(f"| `{study}` | {n} |")
    lines.extend(
        [
            "",
            "### Why the census gap remains",
            "",
            remote.get("note", ""),
            "",
            f"- Remote probe empty-on-both-sides: **{remote.get('n_empty_both_sides')}**",
            f"- Remote recoverable (delta>0): **{remote.get('n_recoverable_delta_gt0')}**",
            "",
            "## Study `platform_id` nulls",
            "",
            f"- Still null: **{plat['n_studies']}** studies (**{plat['n_samples']:,}** samples)",
            "",
            "| Reason | Studies |",
            "| --- | ---: |",
        ]
    )
    for reason, n in plat["reason_counts"].items():
        lines.append(f"| `{reason}` | {n} |")
    lines.extend(
        [
            "",
            "- `mixed_census`: Hub census lists ≥2 Illumina platforms for the study — "
            "study-level `platform_id` correctly stays null (`metadata_json.datahub.platforms` "
            "holds the set).",
            "- `assay_probe_count`: fillable from EWAS_db `*.txt` line counts (HM450/EPIC/EPICv2 bands).",
            "- `no_single_platform_source`: no census/GEO/assay agreement.",
            "",
            "## Notes",
            "",
            "- Planned empty-dir refill lists (`configs/data/ewas_db_refill_*.txt`) are fully "
            "applied; remaining empty dirs have **no** Hub remote `*.txt`.",
            "- Non-GSM namespaces (TCGA/CPTAC/ArrayExpress/…) are mirrored where Hub hosts files.",
            "- Re-run: `uv run python scripts/write_residual_gaps_report.py [--apply-platform-fill]`",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--apply-platform-fill",
        action="store_true",
        help="Update release study.platform_id from EWAS_db probe-count bands",
    )
    parser.add_argument("--release-id", default=RELEASE_ID)
    args = parser.parse_args()
    paths = DataPaths.from_environment()
    report_dir = paths.project_root / "reports" / "inspection" / "deepmat_data_v1"
    report_dir.mkdir(parents=True, exist_ok=True)
    summary = build_residual_gaps(
        paths=paths,
        release_id=args.release_id,
        apply_platform_fill=args.apply_platform_fill,
    )
    write_json(report_dir / "residual_gaps.json", summary)
    write_markdown(summary, report_dir / "residual_gaps.md")
    classification = pd.DataFrame(summary["platform_null"]["classification"])
    classification.to_csv(report_dir / "platform_null_classification.csv", index=False)
    print(
        json.dumps(
            {
                "report_dir": str(report_dir),
                "n_catalog_samples": summary["catalog"]["n_samples"],
                "census_not_in_catalog": summary["vs_datahub_census"]["census_not_in_catalog"],
                "platform_null_studies": summary["platform_null"]["n_studies"],
                "platform_null_samples": summary["platform_null"]["n_samples"],
                "assay_platform_fill": summary["assay_platform_fill"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
