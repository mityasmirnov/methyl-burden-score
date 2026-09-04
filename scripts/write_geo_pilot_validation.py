#!/usr/bin/env python3
"""Write repaired-pilot validation artifacts (tissue/age/conflicts/per-GSE).

Run after rebuilding the 15-GSE parquet with current code and a clean
catalog refresh (zero GEO → merge). Does not expand beyond the pilot list.
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import duckdb
import pandas as pd

from mbs.annotation.manifest import utc_now_iso, write_json
from mbs.geo_metadata import GEO_SOURCE_FAMILY, geo_parquet_path, load_geo_frame
from mbs.paths import DataPaths

PILOT_STUDIES = Path("configs/data/geo_backfill_pilot_gse.txt")
REPORT_DIR = Path("reports/inspection/deepmat_data_v1/geo_backfill_pilot")


def _study_ids(path: Path) -> list[str]:
    return [
        line.strip().upper()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith("#")
    ]


def main() -> None:
    paths = DataPaths.from_environment()
    report_dir = paths.project_root / REPORT_DIR
    report_dir.mkdir(parents=True, exist_ok=True)

    frame = load_geo_frame(paths.data_root)
    if frame.empty:
        raise SystemExit(f"missing GEO parquet: {geo_parquet_path(paths.data_root)}")

    pilot = set(_study_ids(paths.project_root / PILOT_STUDIES))
    studies_in_parquet = set(frame["study_id"].astype(str).str.upper())
    if not studies_in_parquet.issubset(pilot):
        extra = sorted(studies_in_parquet - pilot)
        raise SystemExit(
            f"parquet contains non-pilot studies {extra}; rebuild with pilot list first"
        )

    # --- parquet-side audits (repaired parser) ---
    tissue_by_study: dict[str, dict[str, int]] = {}
    age_by_study: dict[str, dict[str, object]] = {}
    pheno_by_study: dict[str, dict[str, int]] = {}
    for study_id, grp in frame.groupby(frame["study_id"].astype(str), sort=True):
        tmap = Counter(grp.get("tissue_map_status", pd.Series(dtype=str)).fillna("empty"))
        tissue_by_study[str(study_id)] = {
            "mapped": int(tmap.get("mapped", 0)),
            "unmapped": int(tmap.get("unmapped", 0)),
            "ambiguous": int(tmap.get("ambiguous", 0)),
            "empty": int(tmap.get("empty", 0)),
            "n_gsm": int(len(grp)),
            "top_unmapped": [
                {"label": str(lab), "n": int(n)}
                for lab, n in Counter(
                    grp.loc[grp.get("tissue_map_status") == "unmapped", "tissue_raw"]
                    .dropna()
                    .astype(str)
                ).most_common(5)
            ],
            "top_mapped": [
                {"label": str(lab), "n": int(n)}
                for lab, n in Counter(
                    grp.loc[grp.get("tissue_map_status") == "mapped", "tissue"]
                    .dropna()
                    .astype(str)
                ).most_common(5)
            ],
        }
        ages = grp.loc[grp["age"].notna()]
        units = Counter(ages.get("age_unit", pd.Series(dtype=str)).fillna("unknown"))
        age_by_study[str(study_id)] = {
            "n_age": int(len(ages)),
            "age_unit_counts": {str(k): int(v) for k, v in units.items()},
            "age_years_min": float(ages["age"].min()) if len(ages) else None,
            "age_years_median": float(ages["age"].median()) if len(ages) else None,
            "age_years_max": float(ages["age"].max()) if len(ages) else None,
            "age_raw_examples": [
                str(x) for x in ages.get("age_raw", pd.Series(dtype=str)).dropna().astype(str).head(5)
            ],
        }
        pheno_by_study[str(study_id)] = {
            "age": int(grp["age"].notna().sum()),
            "sex": int(grp["sex"].notna().sum()),
            "tissue_mapped": int((grp.get("tissue_map_status") == "mapped").sum()),
            "disease": int(grp["disease"].notna().sum()) if "disease" in grp.columns else 0,
            "cancer": int(grp["cancer"].notna().sum()) if "cancer" in grp.columns else 0,
        }

    multi = 0
    if "study_ids" in frame.columns:
        for raw in frame["study_ids"].tolist():
            try:
                ids = json.loads(raw) if isinstance(raw, str) else raw
            except json.JSONDecodeError:
                ids = []
            if isinstance(ids, list) and len(ids) > 1:
                multi += 1

    fetch_status_path = report_dir / "fetch_status.json"
    fetch_status = (
        json.loads(fetch_status_path.read_text(encoding="utf-8"))
        if fetch_status_path.is_file()
        else None
    )

    catalog = paths.data_root / "canonical" / "releases" / "deepmat-data-v1" / "catalog" / "catalog.duckdb"
    catalog_stats: dict[str, object] = {"catalog_present": catalog.is_file()}
    if catalog.is_file():
        con = duckdb.connect(str(catalog), read_only=True)
        try:
            n_geo = con.execute(
                "SELECT count(*) FROM sample_phenotype WHERE source_family = ?",
                [GEO_SOURCE_FAMILY],
            ).fetchone()
            catalog_stats["n_geo_phenotype_rows"] = int(n_geo[0]) if n_geo else 0
            by_id = con.execute(
                """
                SELECT phenotype_id, label_status, count(*) AS n
                FROM sample_phenotype
                WHERE source_family = ?
                GROUP BY 1, 2
                ORDER BY 1, 2
                """,
                [GEO_SOURCE_FAMILY],
            ).fetchdf()
            catalog_stats["label_status"] = by_id.to_dict(orient="records")
            disease_cases = con.execute(
                """
                SELECT count(*) FROM sample_phenotype
                WHERE source_family = ? AND phenotype_id IN ('disease','cancer')
                  AND label_status = 'case' AND is_observed
                """,
                [GEO_SOURCE_FAMILY],
            ).fetchone()
            catalog_stats["n_disease_cancer_cases"] = int(disease_cases[0]) if disease_cases else 0
            hub_overlap = con.execute(
                """
                SELECT count(*) FROM sample_phenotype sp
                JOIN sample_source_membership m USING (sample_id)
                WHERE sp.source_family = ?
                """,
                [GEO_SOURCE_FAMILY],
            ).fetchone()
            catalog_stats["hub_overlap_geo_rows"] = int(hub_overlap[0]) if hub_overlap else 0
        finally:
            con.close()

    age_units_global = Counter(
        frame.loc[frame["age"].notna(), "age_unit"].fillna("unknown").astype(str)
    )
    tissue_global = Counter(frame.get("tissue_map_status", pd.Series(dtype=str)).fillna("empty"))

    summary = {
        "generated_at": utc_now_iso(),
        "validation_kind": "repaired_pilot_15gse",
        "n_parquet_gsm": int(len(frame)),
        "study_ids": sorted(studies_in_parquet),
        "n_multi_study_gsm": multi,
        "conflict_stats_from_fetch": (fetch_status or {}).get("conflict_stats"),
        "conflicts_from_fetch": (fetch_status or {}).get("conflicts") or [],
        "tissue_map_global": {str(k): int(v) for k, v in tissue_global.items()},
        "age_unit_global": {str(k): int(v) for k, v in age_units_global.items()},
        "age_years_global": {
            "n": int(frame["age"].notna().sum()),
            "min": float(frame["age"].min()) if frame["age"].notna().any() else None,
            "median": float(frame["age"].median()) if frame["age"].notna().any() else None,
            "max": float(frame["age"].max()) if frame["age"].notna().any() else None,
        },
        "tissue_by_study": tissue_by_study,
        "age_by_study": age_by_study,
        "phenotype_counts_by_study": pheno_by_study,
        "catalog": catalog_stats,
        "training_notes": [
            "Training heads read Hub pack Parquet, not geo_metadata_backfill.",
            "GEO improves catalog / eligibility only until a separate geo-dev release.",
            "Disease/cancer cases must be >0 with eligibility cutoffs before any disease head.",
            "Do not mutate frozen ATS (matrix-hub-age-tissue-sex-full-v1).",
        ],
    }
    write_json(report_dir / "validation.json", summary)

    lines = [
        "# GEO repaired-pilot validation (15 GSE)",
        "",
        f"- Generated: `{summary['generated_at']}`",
        f"- Parquet GSM: **{summary['n_parquet_gsm']}**",
        f"- Multi-study GSM (membership): **{multi}**",
        f"- Fetch conflict samples: **{(summary.get('conflict_stats_from_fetch') or {}).get('n_conflict_samples', 'n/a')}**",
        "",
        "## Global tissue map",
        "",
        "| status | n |",
        "| --- | ---: |",
    ]
    for k, v in sorted(summary["tissue_map_global"].items()):
        lines.append(f"| `{k}` | {v} |")
    lines.extend(
        [
            "",
            "## Global age units (after conversion to years)",
            "",
            "| age_unit | n |",
            "| --- | ---: |",
        ]
    )
    for k, v in sorted(summary["age_unit_global"].items()):
        lines.append(f"| `{k}` | {v} |")
    ag = summary["age_years_global"]
    lines.extend(
        [
            "",
            f"- Age years min/median/max: **{ag['min']}** / **{ag['median']}** / **{ag['max']}** (n={ag['n']})",
            "",
            "## Catalog (after merge)",
            "",
            f"- GEO phenotype rows: **{catalog_stats.get('n_geo_phenotype_rows', 'n/a')}**",
            f"- Disease/cancer **cases**: **{catalog_stats.get('n_disease_cancer_cases', 'n/a')}** (must not train disease head if 0)",
            f"- GEO rows on Hub membership GSM: **{catalog_stats.get('hub_overlap_geo_rows', 'n/a')}** (must be 0)",
            "",
            "## Tissue coverage by study",
            "",
            "| study_id | gsm | mapped | unmapped | empty | top unmapped |",
            "| --- | ---: | ---: | ---: | ---: | --- |",
        ]
    )
    for sid, rec in tissue_by_study.items():
        top_u = ", ".join(f"{x['label']} ({x['n']})" for x in rec["top_unmapped"][:3]) or "—"
        lines.append(
            f"| `{sid}` | {rec['n_gsm']} | {rec['mapped']} | {rec['unmapped']} | "
            f"{rec['empty']} | {top_u} |"
        )
    lines.extend(
        [
            "",
            "## Operator notes",
            "",
        ]
    )
    lines.extend(f"- {n}" for n in summary["training_notes"])
    lines.append("")
    (report_dir / "validation.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {report_dir / 'validation.json'}")
    print(f"wrote {report_dir / 'validation.md'}")


if __name__ == "__main__":
    main()
