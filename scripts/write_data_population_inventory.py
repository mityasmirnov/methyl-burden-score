#!/usr/bin/env python3
"""Write deepmat-data-v1 data population inventory (catalog + EWAS_db disk + GEO)."""

from __future__ import annotations

import re
from pathlib import Path

import duckdb
import pandas as pd

from mbs.annotation.manifest import utc_now_iso, write_json
from mbs.paths import DataPaths


def main() -> None:
    paths = DataPaths.from_environment()
    con = duckdb.connect(
        str(paths.data_root / "canonical/releases/deepmat-data-v1/catalog/catalog.duckdb"),
        read_only=True,
    )

    def q(sql: str) -> pd.DataFrame:
        return con.execute(sql).fetchdf()

    out: dict = {"generated_at": utc_now_iso(), "release_id": "deepmat-data-v1"}
    out["catalog"] = {
        "n_samples": int(q("SELECT count(*) n FROM sample")["n"][0]),
        "n_studies": int(q("SELECT count(*) n FROM study")["n"][0]),
        "n_phenotype_rows": int(q("SELECT count(*) n FROM sample_phenotype")["n"][0]),
        "n_assay_files": int(q("SELECT count(*) n FROM assay_file")["n"][0]),
    }
    out["platforms"] = q(
        """
        SELECT coalesce(s.platform_id, '(null)') AS platform_id, count(sa.sample_id) AS n_samples
        FROM sample sa LEFT JOIN study s ON sa.study_id = s.study_id
        GROUP BY 1 ORDER BY 2 DESC
        """
    ).to_dict(orient="records")
    out["phenotypes"] = q(
        """
        SELECT phenotype_id, count(*) n_rows, count(DISTINCT sample_id) n_gsm,
          sum(CASE WHEN numeric_value IS NOT NULL THEN 1 ELSE 0 END) n_numeric,
          sum(CASE WHEN categorical_value IS NOT NULL AND categorical_value != '' THEN 1 ELSE 0 END)
            n_categorical
        FROM sample_phenotype GROUP BY 1 ORDER BY n_gsm DESC
        """
    ).to_dict(orient="records")
    out["by_source_family"] = q(
        """
        SELECT source_family, phenotype_id, count(DISTINCT sample_id) n_gsm, count(*) n_rows
        FROM sample_phenotype GROUP BY 1,2 ORDER BY n_gsm DESC
        """
    ).to_dict(orient="records")
    out["tissue_top"] = q(
        """
        SELECT categorical_value AS tissue, count(DISTINCT sample_id) n_gsm
        FROM sample_phenotype
        WHERE phenotype_id='tissue' AND categorical_value IS NOT NULL AND categorical_value != ''
        GROUP BY 1 ORDER BY 2 DESC LIMIT 40
        """
    ).to_dict(orient="records")
    out["age_sex"] = {
        "age_gsm": int(
            q(
                "SELECT count(DISTINCT sample_id) n FROM sample_phenotype "
                "WHERE phenotype_id='age' AND numeric_value IS NOT NULL"
            )["n"][0]
        ),
        "sex_gsm": int(
            q("SELECT count(DISTINCT sample_id) n FROM sample_phenotype WHERE phenotype_id='sex'")[
                "n"
            ][0]
        ),
        "tissue_gsm": int(
            q(
                "SELECT count(DISTINCT sample_id) n FROM sample_phenotype "
                "WHERE phenotype_id='tissue' AND categorical_value IS NOT NULL"
            )["n"][0]
        ),
    }
    for pid in ("disease", "cancer"):
        out[f"{pid}_label_status"] = q(
            f"""
            SELECT coalesce(label_status,'(null)') label_status, count(DISTINCT sample_id) n_gsm
            FROM sample_phenotype WHERE phenotype_id='{pid}' GROUP BY 1 ORDER BY 2 DESC
            """
        ).to_dict(orient="records")
    out["trait_eligibility"] = q(
        """
        SELECT phenotype_id, phenotype_family, task_type, n_samples, n_cases, n_controls, n_unknown,
               prevalence, n_studies, n_platforms, n_tissues,
               eligible_core_task, eligible_auxiliary_task, eligible_external_evaluation,
               exclusion_reason
        FROM trait_eligibility ORDER BY phenotype_family, phenotype_id
        """
    ).to_dict(orient="records")
    out["geo_series"] = q(
        """
        SELECT
          count(*) FILTER (WHERE study_id LIKE 'GSE%') AS n_gse,
          count(*) FILTER (
            WHERE json_extract_string(metadata_json,'$.geo.title') IS NOT NULL
          ) AS with_title,
          count(*) FILTER (
            WHERE length(trim(coalesce(
              json_extract_string(metadata_json,'$.geo.overall_design'),'')))>0
          ) AS with_overall_design
        FROM study
        """
    ).to_dict(orient="records")[0]
    con.close()

    ewas = paths.data_root / "raw" / "ewas_datahub" / "EWAS_db"
    gsm_re = re.compile(r"^GSM[0-9]+\.txt$", re.I)
    txt_re = re.compile(r"\.txt$", re.I)
    n_dirs = n_gsm_dirs = n_txt_dirs = gsm = txt = 0
    for d in ewas.iterdir():
        if not d.is_dir():
            continue
        n_dirs += 1
        files = [p.name for p in d.iterdir() if p.is_file()]
        ng = sum(1 for n in files if gsm_re.match(n))
        nt = sum(1 for n in files if txt_re.search(n))
        gsm += ng
        txt += nt
        if ng:
            n_gsm_dirs += 1
        if nt:
            n_txt_dirs += 1
    out["ewas_db_disk"] = {
        "n_study_dirs": n_dirs,
        "n_dirs_with_gsm_txt": n_gsm_dirs,
        "n_dirs_with_any_txt": n_txt_dirs,
        "n_gsm_txt_files": gsm,
        "n_txt_files": txt,
        "advertised_studies": 1989,
    }
    gdf = pd.read_parquet(
        paths.data_root / "canonical/phenotypes/geo_sample_metadata.parquet",
        columns=["study_id", "age", "sex", "tissue_map_status"],
    )
    out["geo_parquet"] = {
        "n_gsm": len(gdf),
        "n_studies": int(gdf["study_id"].nunique()),
        "age": int(gdf["age"].notna().sum()),
        "sex": int(gdf["sex"].notna().sum()),
        "tissue_map": {
            str(k): int(v) for k, v in gdf["tissue_map_status"].value_counts(dropna=False).items()
        },
    }

    report_dir = Path("reports/inspection/deepmat_data_v1")
    write_json(report_dir / "data_population_inventory.json", out)

    def yn(b: object) -> str:
        return "yes" if b else "no"

    lines = [
        "# Data population inventory",
        "",
        f"- Generated: `{out['generated_at']}`",
        f"- Release: `{out['release_id']}`",
        "",
        "## Catalog scale",
        "",
        "| Samples | Studies | Phenotype rows | Assay files |",
        "| ---: | ---: | ---: | ---: |",
        f"| **{out['catalog']['n_samples']:,}** | **{out['catalog']['n_studies']:,}** | "
        f"**{out['catalog']['n_phenotype_rows']:,}** | **{out['catalog']['n_assay_files']:,}** |",
        "",
        "## Assays / platforms",
        "",
        "Platform is study-level (`study.platform_id`); counts are samples under those studies.",
        "",
        "| Platform | Samples |",
        "| --- | ---: |",
    ]
    for r in out["platforms"]:
        lines.append(f"| `{r['platform_id']}` | {int(r['n_samples']):,} |")
    lines += [
        "",
        "## Phenotypes / traits",
        "",
        "| phenotype_id | distinct GSM | rows | numeric | categorical |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for r in out["phenotypes"]:
        lines.append(
            f"| `{r['phenotype_id']}` | {int(r['n_gsm']):,} | {int(r['n_rows']):,} | "
            f"{int(r['n_numeric']):,} | {int(r['n_categorical']):,} |"
        )
    lines += [
        "",
        f"- Age (numeric): **{out['age_sex']['age_gsm']:,}** GSM",
        f"- Sex: **{out['age_sex']['sex_gsm']:,}** GSM",
        f"- Tissue (mapped categorical in catalog): **{out['age_sex']['tissue_gsm']:,}** GSM",
        "",
        "### Top tissues",
        "",
        "| Tissue | GSM |",
        "| --- | ---: |",
    ]
    for r in out["tissue_top"][:25]:
        lines.append(f"| {r['tissue']} | {int(r['n_gsm']):,} |")
    lines += ["", "### Disease / cancer label status", ""]
    for pid in ("disease", "cancer"):
        lines += [f"**{pid}**", "", "| label_status | GSM |", "| --- | ---: |"]
        for r in out[f"{pid}_label_status"]:
            lines.append(f"| `{r['label_status']}` | {int(r['n_gsm']):,} |")
        lines.append("")
    lines += [
        "## Trait eligibility (training gates)",
        "",
        "| phenotype | family | task | n | cases | controls | unknown | core | aux | external | exclusion |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | --- | --- | --- | --- |",
    ]
    for r in out["trait_eligibility"]:
        lines.append(
            f"| `{r['phenotype_id']}` | {r['phenotype_family']} | {r['task_type']} | "
            f"{r['n_samples']} | {r['n_cases']} | {r['n_controls']} | {r['n_unknown']} | "
            f"{yn(r['eligible_core_task'])} | {yn(r['eligible_auxiliary_task'])} | "
            f"{yn(r['eligible_external_evaluation'])} | {r['exclusion_reason'] or ''} |"
        )
    lines += [
        "",
        "## Phenotype provenance (top source_family × phenotype)",
        "",
        "| source_family | phenotype_id | GSM | rows |",
        "| --- | --- | ---: | ---: |",
    ]
    for r in out["by_source_family"][:40]:
        lines.append(
            f"| `{r['source_family']}` | `{r['phenotype_id']}` | "
            f"{int(r['n_gsm']):,} | {int(r['n_rows']):,} |"
        )
    ed, g, gs = out["ewas_db_disk"], out["geo_parquet"], out["geo_series"]
    lines += [
        "",
        "## EWAS_db assay mirror (on-disk methylation profiles)",
        "",
        "| Advertised | Local dirs | Dirs with GSM*.txt | GSM*.txt | Any *.txt |",
        "| ---: | ---: | ---: | ---: | ---: |",
        f"| {ed['advertised_studies']} | {ed['n_study_dirs']} | {ed['n_dirs_with_gsm_txt']} | "
        f"{ed['n_gsm_txt_files']:,} | {ed['n_txt_files']:,} |",
        "",
        "Empty-dir refill: `scripts/refill_ewas_db_empty_studies.sh` + "
        "`configs/data/ewas_db_refill_*.txt` (see `ewas_db_empty_refill_plan.json`).",
        "",
        "## GEO metadata (labels / study context — not assay betas)",
        "",
        f"- Sample parquet: **{g['n_gsm']:,}** GSM / **{g['n_studies']}** primary studies; "
        f"age **{g['age']:,}**, sex **{g['sex']:,}**, tissue_map `{g['tissue_map']}`",
        f"- Series: **{gs['n_gse']}** GSE; title **{gs['with_title']}**; "
        f"overall_design **{gs['with_overall_design']}**",
        "",
        "## Related reports",
        "",
        "- [`census.md`](census.md)",
        "- [`trait_eligibility.md`](trait_eligibility.md)",
        "- [`geo_series_enrichment.md`](geo_series_enrichment.md)",
        "- [`ewas_db_download_failures.md`](ewas_db_download_failures.md)",
        "",
    ]
    (report_dir / "data_population_inventory.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {report_dir / 'data_population_inventory.md'}")


if __name__ == "__main__":
    main()
