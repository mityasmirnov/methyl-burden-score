#!/usr/bin/env python3
"""Write sample overview inspection report, figures, and PDF dashboard."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, cast

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

from mbs.annotation.manifest import utc_now_iso, write_json
from mbs.paths import DataPaths
from mbs.sample_overview_enrich import CATALOG_PLATFORMS, platform_from_n_probes

HUB_FAMILIES = (
    "age",
    "tissue",
    "sex",
    "disease",
    "cancer",
    "blood",
    "brain",
    "bmi",
    "ancestry",
)

_DISEASE_DISPLAY = {
    "alzheimer's disease": "Alzheimer's disease",
    "parkinson's disease": "Parkinson's disease",
    "huntington's disease": "Huntington's disease",
    "crohn's disease": "Crohn's disease",
    "graves' disease": "Graves' disease",
    "sjogren's syndrome": "Sjogren's syndrome",
    "down syndrome": "Down syndrome",
    "silver russell syndrome": "Silver Russell syndrome",
    "kabuki syndrome": "Kabuki syndrome",
}


def _norm_disease(raw: object) -> str:
    text = str(raw).strip().lower().replace("\u2019", "'")
    if text in {"crohns disease", "crohn's disease"}:
        return "crohn's disease"
    if text == "ulcerative colitis":
        return "ulcerative colitis"
    return text


def _pretty_disease(key: str) -> str:
    return _DISEASE_DISPLAY.get(key, key)


def build_disease_census(*, data_root: Path) -> dict[str, object]:
    """Hub disease-pack census for the nine-pack cohort (case/control + top diagnoses)."""
    phenotypes = data_root / "canonical" / "phenotypes"
    label_path = phenotypes / "sample_phenotype_table_hub_nine_pack_v1.label_status.parquet"
    info_path = phenotypes / "disease_sample_info.parquet"
    if not label_path.is_file() or not info_path.is_file():
        return {}

    labels = pd.read_parquet(
        label_path,
        columns=[
            "sample_id",
            "disease_mask",
            "disease_label_status",
            "disease_is_case",
            "disease_case_control_mask",
            "cancer_mask",
            "cancer_label_status",
            "cancer_is_case",
        ],
    )
    n_nine = len(labels)
    pack = labels.loc[labels["disease_mask"]].drop_duplicates("sample_id")
    n_pack = len(pack)
    n_case = int(pack["disease_is_case"].sum())
    n_control = int((pack["disease_label_status"] == "control").sum())
    n_adjacent = int((pack["disease_label_status"] == "adjacent_normal").sum())
    n_usable = int(pack["disease_case_control_mask"].sum())

    info = pd.read_parquet(info_path, columns=["sample_id", "phenotype_value"])
    cases = info.merge(pack[["sample_id", "disease_is_case"]], on="sample_id", how="inner")
    cases = cases.loc[cases["disease_is_case"]].drop_duplicates("sample_id")
    cases = cases.copy()
    cases["disease_key"] = cases["phenotype_value"].map(_norm_disease)
    counts = cases["disease_key"].value_counts()
    top: list[dict[str, object]] = []
    for key, n in counts.items():
        top.append(
            {
                "disease": _pretty_disease(str(key)),
                "n_cases": int(n),
                "pct_of_cases": round(100.0 * float(n) / float(n_case), 2) if n_case else 0.0,
            }
        )

    cancer_pack = labels.loc[labels["cancer_mask"]].drop_duplicates("sample_id")
    cancer = {
        "n_pack": int(len(cancer_pack)),
        "n_case": int(cancer_pack["cancer_is_case"].sum()),
        "n_control": int((cancer_pack["cancer_label_status"] == "control").sum()),
        "n_adjacent_normal": int((cancer_pack["cancer_label_status"] == "adjacent_normal").sum()),
        "frac_patients_of_pack": (
            round(float(cancer_pack["cancer_is_case"].mean()), 4) if len(cancer_pack) else 0.0
        ),
    }

    return {
        "n_nine_pack": n_nine,
        "n_pack": n_pack,
        "n_case": n_case,
        "n_control": n_control,
        "n_adjacent_normal": n_adjacent,
        "n_usable_case_control": n_usable,
        "n_disease_labels": len(top),
        "frac_patients_of_pack": round(n_case / n_pack, 4) if n_pack else 0.0,
        "frac_patients_of_usable": round(n_case / n_usable, 4) if n_usable else 0.0,
        "frac_patients_of_ninepack": round(n_case / n_nine, 4) if n_nine else 0.0,
        "top_diseases": top,
        "cancer": cancer,
        "note": (
            "Disease names come from Hub disease-pack phenotype_value on case rows; "
            "controls have sample_type=control and usually no diagnosis string. "
            "Train only on label_status in {case, control}."
        ),
    }


def _coverage(series: pd.Series) -> int:
    return int(series.notna().sum())


def _as_dict(value: object, *, name: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise TypeError(f"expected dict for {name}, got {type(value).__name__}")
    return cast(dict[str, Any], value)


def build_summary(frame: pd.DataFrame, *, data_root: Path) -> dict[str, object]:
    n = len(frame)
    hub = int(frame["in_hub_baseline"].fillna(False).sum()) if "in_hub_baseline" in frame else 0
    ewas = int(frame["in_ewas_db"].fillna(False).sum()) if "in_ewas_db" in frame else 0
    both = int(
        ((frame["in_hub_baseline"].fillna(False)) & (frame["in_ewas_db"].fillna(False))).sum()
    )
    geo = int(frame["in_geo_soft"].fillna(False).sum()) if "in_geo_soft" in frame else 0
    pack_counts = {
        f: int(frame[f"in_hub_{f}"].fillna(False).sum())
        for f in HUB_FAMILIES
        if f"in_hub_{f}" in frame.columns
    }
    platform = (
        frame["platform_id"].fillna("unknown").astype(str).value_counts().head(12).to_dict()
        if "platform_id" in frame.columns
        else {}
    )
    platform_source = (
        frame["platform_source"].fillna("missing").astype(str).value_counts().to_dict()
        if "platform_source" in frame.columns
        else {}
    )
    species = (
        frame["species_status"].fillna("missing").astype(str).value_counts().to_dict()
        if "species_status" in frame.columns
        else {}
    )
    species_source = (
        frame["species_source"].fillna("missing").astype(str).value_counts().to_dict()
        if "species_source" in frame.columns
        else {}
    )
    phenotype_coverage = {
        "age_years": _coverage(frame["age_years"]) if "age_years" in frame else 0,
        "sex": _coverage(frame["sex"]) if "sex" in frame else 0,
        "tissue": _coverage(frame["tissue"]) if "tissue" in frame else 0,
        "bmi": _coverage(frame["bmi"]) if "bmi" in frame else 0,
        "ancestry_label": _coverage(frame["ancestry_label"]) if "ancestry_label" in frame else 0,
        "brain_label": _coverage(frame["brain_label"]) if "brain_label" in frame else 0,
        "cell_component": _coverage(frame["cell_component"]) if "cell_component" in frame else 0,
        "pubmed_ids": _coverage(frame["pubmed_ids"]) if "pubmed_ids" in frame else 0,
    }
    sex_top = (
        frame["sex"].dropna().astype(str).value_counts().head(10).to_dict()
        if "sex" in frame.columns
        else {}
    )
    top_studies = (
        frame["study_id"].astype(str).value_counts().head(15).to_dict()
        if "study_id" in frame.columns
        else {}
    )
    tissue_top = (
        frame["tissue"].dropna().astype(str).value_counts().head(20).to_dict()
        if "tissue" in frame.columns
        else {}
    )
    ancestry_top = (
        frame["ancestry_label"].dropna().astype(str).value_counts().head(15).to_dict()
        if "ancestry_label" in frame.columns
        else {}
    )
    brain_top = (
        frame["brain_label"].dropna().astype(str).value_counts().head(12).to_dict()
        if "brain_label" in frame.columns
        else {}
    )
    disease = build_disease_census(data_root=data_root)
    return {
        "created_at": utc_now_iso(),
        "n_samples": n,
        "lane": {
            "in_hub_baseline": hub,
            "in_ewas_db": ewas,
            "hub_and_ewas": both,
            "in_geo_soft": geo,
        },
        "hub_pack_membership": pack_counts,
        "platform_top": platform,
        "platform_source": platform_source,
        "species_status": species,
        "species_source": species_source,
        "phenotype_coverage": phenotype_coverage,
        "sex_top": sex_top,
        "top_studies": top_studies,
        "tissue_top": tissue_top,
        "ancestry_top": ancestry_top,
        "brain_top": brain_top,
        "disease": disease,
        "cpgpt_sample_embedding_status": "unavailable",
        "id_convention": {
            "sample_id": "GSM (or TCGA/ArrayExpress when not GEO)",
            "study_id": "GSE (series)",
            "geo_gpl_id": "GPL (platform)",
        },
    }


def _barh(ax: Any, mapping: dict[str, Any], *, title: str, xlabel: str) -> None:
    keys = list(mapping.keys())
    vals = [int(mapping[k]) for k in keys]
    ax.barh(keys[::-1], vals[::-1], color="#4C78A8")
    ax.set_xlabel(xlabel)
    ax.set_title(title)


def write_figures(summary: dict[str, object], fig_dir: Path) -> list[str]:
    fig_dir.mkdir(parents=True, exist_ok=True)
    written: list[str] = []

    lane = _as_dict(summary["lane"], name="lane")
    fig, ax = plt.subplots(figsize=(7, 4))
    labels = ["Hub baseline", "EWAS_db", "Hub and EWAS", "GEO SOFT"]
    values = [
        int(lane["in_hub_baseline"]),
        int(lane["in_ewas_db"]),
        int(lane["hub_and_ewas"]),
        int(lane["in_geo_soft"]),
    ]
    ax.bar(labels, values, color="#4C78A8")
    ax.set_ylabel("Samples")
    ax.set_title("Catalog lane membership")
    ax.tick_params(axis="x", rotation=20)
    fig.tight_layout()
    path = fig_dir / "lane_membership.png"
    fig.savefig(path, dpi=120)
    plt.close(fig)
    written.append(str(path))

    packs = _as_dict(summary["hub_pack_membership"], name="hub_pack_membership")
    fig, ax = plt.subplots(figsize=(8, 4))
    keys = list(packs.keys())
    vals = [int(packs[k]) for k in keys]
    ax.bar(keys, vals, color="#F58518")
    ax.set_ylabel("Samples")
    ax.set_title("Hub pack membership (catalog)")
    ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    path = fig_dir / "hub_pack_membership.png"
    fig.savefig(path, dpi=120)
    plt.close(fig)
    written.append(str(path))

    cov = _as_dict(summary["phenotype_coverage"], name="phenotype_coverage")
    fig, ax = plt.subplots(figsize=(8, 4))
    _barh(ax, cov, title="Phenotype / metadata coverage", xlabel="Samples with non-null value")
    fig.tight_layout()
    path = fig_dir / "phenotype_coverage.png"
    fig.savefig(path, dpi=120)
    plt.close(fig)
    written.append(str(path))

    tissue = _as_dict(summary["tissue_top"], name="tissue_top")
    fig, ax = plt.subplots(figsize=(8, 6))
    _barh(ax, tissue, title="Top tissues (harmonized label)", xlabel="Samples")
    fig.tight_layout()
    path = fig_dir / "tissue_top.png"
    fig.savefig(path, dpi=120)
    plt.close(fig)
    written.append(str(path))

    platform = _as_dict(summary["platform_top"], name="platform_top")
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(list(platform.keys()), [int(v) for v in platform.values()], color="#54A24B")
    ax.set_ylabel("Samples")
    ax.set_title("Platform (after inference)")
    fig.tight_layout()
    path = fig_dir / "platform.png"
    fig.savefig(path, dpi=120)
    plt.close(fig)
    written.append(str(path))

    disease = summary.get("disease")
    if isinstance(disease, dict) and disease.get("top_diseases"):
        top = disease["top_diseases"]
        assert isinstance(top, list)
        mapping = {
            str(row["disease"]): int(row["n_cases"])  # type: ignore[index]
            for row in top[:15]
            if isinstance(row, dict)
        }
        fig, ax = plt.subplots(figsize=(8, 6))
        _barh(ax, mapping, title="Top disease diagnoses (cases)", xlabel="Case samples")
        fig.tight_layout()
        path = fig_dir / "disease_top.png"
        fig.savefig(path, dpi=120)
        plt.close(fig)
        written.append(str(path))

    return written


def write_pdf_dashboard(summary: dict[str, object], path: Path) -> None:
    """Multi-page PDF suitable for GitHub (reports/inspection/...)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lane = _as_dict(summary["lane"], name="lane")
    packs = _as_dict(summary["hub_pack_membership"], name="hub_pack_membership")
    cov = _as_dict(summary["phenotype_coverage"], name="phenotype_coverage")
    tissue = _as_dict(summary["tissue_top"], name="tissue_top")
    platform = _as_dict(summary["platform_top"], name="platform_top")
    species = _as_dict(summary["species_status"], name="species_status")
    species_source = _as_dict(summary["species_source"], name="species_source")
    platform_source = _as_dict(summary["platform_source"], name="platform_source")
    sex_top = _as_dict(summary.get("sex_top", {}), name="sex_top")
    ancestry = _as_dict(summary.get("ancestry_top", {}), name="ancestry_top")
    brain = _as_dict(summary.get("brain_top", {}), name="brain_top")

    with PdfPages(path) as pdf:
        # Page 1 — headline
        fig = plt.figure(figsize=(11, 8.5))
        fig.suptitle("Sample overview — Hub + EWAS_db", fontsize=16, fontweight="bold")
        fig.text(
            0.5,
            0.92,
            f"n={summary['n_samples']}  |  {summary['created_at']}  |  "
            "GSM=sample, GSE=study, GPL=platform",
            ha="center",
            fontsize=10,
        )
        ax1 = fig.add_axes([0.08, 0.48, 0.4, 0.35])
        ax1.bar(
            ["Hub", "EWAS_db", "Hub&EWAS", "GEO"],
            [
                int(lane["in_hub_baseline"]),
                int(lane["in_ewas_db"]),
                int(lane["hub_and_ewas"]),
                int(lane["in_geo_soft"]),
            ],
            color="#4C78A8",
        )
        ax1.set_title("Lane membership")
        ax1.set_ylabel("Samples")

        ax2 = fig.add_axes([0.55, 0.48, 0.4, 0.35])
        ax2.bar(list(platform.keys()), [int(v) for v in platform.values()], color="#54A24B")
        ax2.set_title("Platform (inferred)")
        ax2.tick_params(axis="x", rotation=20)

        ax3 = fig.add_axes([0.08, 0.08, 0.4, 0.32])
        ax3.bar(list(species.keys()), [int(v) for v in species.values()], color="#E45756")
        ax3.set_title("Species (all rows filled)")

        ax4 = fig.add_axes([0.55, 0.08, 0.4, 0.32])
        ax4.bar(list(packs.keys()), [int(v) for v in packs.values()], color="#F58518")
        ax4.set_title("Hub pack membership")
        ax4.tick_params(axis="x", rotation=30)
        pdf.savefig(fig)
        plt.close(fig)

        # Page 2 — phenotypes
        fig = plt.figure(figsize=(11, 8.5))
        fig.suptitle("Phenotypes overview", fontsize=16, fontweight="bold")
        ax1 = fig.add_axes([0.12, 0.55, 0.75, 0.35])
        _barh(ax1, cov, title="Coverage (non-null)", xlabel="Samples")
        ax2 = fig.add_axes([0.12, 0.08, 0.35, 0.38])
        if sex_top:
            ax2.bar(list(sex_top.keys()), [int(v) for v in sex_top.values()], color="#72B7B2")
            ax2.set_title("Sex")
            ax2.tick_params(axis="x", rotation=20)
        ax3 = fig.add_axes([0.55, 0.08, 0.4, 0.38])
        if ancestry:
            _barh(ax3, ancestry, title="Ancestry labels", xlabel="Samples")
        pdf.savefig(fig)
        plt.close(fig)

        # Page 3 — tissues
        fig = plt.figure(figsize=(11, 8.5))
        fig.suptitle("Tissues overview", fontsize=16, fontweight="bold")
        ax1 = fig.add_axes([0.2, 0.08, 0.7, 0.8])
        _barh(ax1, tissue, title="Top 20 tissues", xlabel="Samples")
        pdf.savefig(fig)
        plt.close(fig)

        # Page 4 — brain + provenance
        fig = plt.figure(figsize=(11, 8.5))
        fig.suptitle("Brain regions + provenance", fontsize=16, fontweight="bold")
        ax1 = fig.add_axes([0.2, 0.45, 0.7, 0.42])
        if brain:
            _barh(ax1, brain, title="Brain region labels (Hub)", xlabel="Samples")
        else:
            ax1.axis("off")
            ax1.text(0.5, 0.5, "No brain labels", ha="center")
        ax2 = fig.add_axes([0.08, 0.08, 0.4, 0.28])
        ax2.barh(
            list(species_source.keys())[::-1],
            [int(v) for v in list(species_source.values())][::-1],
            color="#B279A2",
        )
        ax2.set_title("Species source")
        ax3 = fig.add_axes([0.55, 0.08, 0.4, 0.28])
        ax3.barh(
            list(platform_source.keys())[::-1],
            [int(v) for v in list(platform_source.values())][::-1],
            color="#FF9D98",
        )
        ax3.set_title("Platform source")
        pdf.savefig(fig)
        plt.close(fig)

        # Page 5 — disease
        disease = summary.get("disease")
        if isinstance(disease, dict) and disease.get("top_diseases"):
            top = cast(list[dict[str, object]], disease["top_diseases"])
            fig = plt.figure(figsize=(11, 8.5))
            fig.suptitle("Disease pack overview (nine-pack Hub)", fontsize=16, fontweight="bold")
            fig.text(
                0.5,
                0.92,
                (
                    f"pack={disease['n_pack']}  cases={disease['n_case']}  "
                    f"controls={disease['n_control']}  "
                    f"patient fraction={float(disease['frac_patients_of_pack']):.1%} of pack  "
                    f"({float(disease['frac_patients_of_ninepack']):.1%} of nine-pack)"
                ),
                ha="center",
                fontsize=10,
            )
            mapping = {
                str(row["disease"]): int(row["n_cases"]) for row in top[:15] if "disease" in row
            }
            ax1 = fig.add_axes([0.22, 0.08, 0.7, 0.78])
            _barh(ax1, mapping, title="Most frequent diagnoses (cases only)", xlabel="Case samples")
            pdf.savefig(fig)
            plt.close(fig)


def write_markdown(summary: dict[str, object], path: Path) -> None:
    lane = _as_dict(summary["lane"], name="lane")
    cov = _as_dict(summary["phenotype_coverage"], name="phenotype_coverage")
    packs = _as_dict(summary["hub_pack_membership"], name="hub_pack_membership")
    platform = _as_dict(summary["platform_top"], name="platform_top")
    species = _as_dict(summary["species_status"], name="species_status")
    lines = [
        "# Sample overview (Hub + EWAS_db)",
        "",
        f"Generated: `{summary['created_at']}`",
        "",
        f"Rows: **{summary['n_samples']}** (one per catalog `sample_id`).",
        "",
        "## ID convention",
        "",
        "- `sample_id` = **GSM** (sample)",
        "- `study_id` / `gse_id` = **GSE** (series/study)",
        "- `geo_gpl_id` = **GPL** (GEO platform)",
        "",
        "## Lane membership",
        "",
        f"- Hub baseline: {lane['in_hub_baseline']}",
        f"- EWAS_db: {lane['in_ewas_db']}",
        f"- Hub and EWAS_db: {lane['hub_and_ewas']}",
        f"- GEO SOFT fetched: {lane['in_geo_soft']}",
        "",
        "## Platform (after inference)",
        "",
        "| Platform | N |",
        "|----------|--:|",
    ]
    for k, v in platform.items():
        lines.append(f"| {k} | {v} |")
    lines.extend(
        [
            "",
            "## Species (all rows filled)",
            "",
            "| Status | N |",
            "|--------|--:|",
        ]
    )
    for k, v in species.items():
        lines.append(f"| {k} | {v} |")
    lines.extend(
        [
            "",
            "## Hub pack membership",
            "",
            "| Pack | N |",
            "|------|--:|",
        ]
    )
    for k, v in packs.items():
        lines.append(f"| {k} | {v} |")
    lines.extend(
        [
            "",
            "## Phenotype coverage (non-null)",
            "",
            "| Field | N |",
            "|-------|--:|",
        ]
    )
    for k, v in cov.items():
        lines.append(f"| {k} | {v} |")
    disease = summary.get("disease")
    if isinstance(disease, dict) and disease:
        lines.extend(
            [
                "",
                "## Disease pack (nine-pack Hub)",
                "",
                f"- Pack membership: **{disease.get('n_pack')}**",
                f"- Cases (patients): **{disease.get('n_case')}** "
                f"({100 * float(disease.get('frac_patients_of_pack', 0)):.1f}% of pack; "
                f"{100 * float(disease.get('frac_patients_of_ninepack', 0)):.1f}% of nine-pack)",
                f"- Controls: **{disease.get('n_control')}**",
                f"- Adjacent normal: **{disease.get('n_adjacent_normal')}**",
                f"- Distinct diagnosis labels (cases): **{disease.get('n_disease_labels')}**",
                "",
                "| Disease | Cases | % of cases |",
                "|---------|------:|----------:|",
            ]
        )
        top = disease.get("top_diseases")
        if isinstance(top, list):
            for row in top:
                if not isinstance(row, dict):
                    continue
                lines.append(
                    f"| {row.get('disease')} | {row.get('n_cases')} | {row.get('pct_of_cases')} |"
                )
        note = disease.get("note")
        if note:
            lines.extend(["", str(note), ""])
    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- Hub phenotype labels win; GEO fills blanks.",
            "- Unknown platforms inferred from DataHub metadata, GPL maps, "
            "`_935k` hints, and EWAS_db probe counts.",
            "- Species filled for every row (GEO SOFT > series taxon > Hub/project "
            "priors > EWAS assumed human); provenance in `species_source`.",
            "- `cpgpt_sample_embedding_status=unavailable` (locus embeddings only).",
            "- Artifacts: `data/canonical/phenotypes/sample_overview_hub_geo_v1.{parquet,csv.gz}`",
            "- PDF dashboard: `sample_overview_dashboard.pdf`",
            "- Canvas (project copy): `sample-overview-hub-geo.canvas.tsx`",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--overview",
        type=Path,
        default=None,
        help="Path to sample_overview_hub_geo_v1.parquet",
    )
    parser.add_argument(
        "--report-dir",
        type=Path,
        default=Path("reports/inspection/deepmat_data_v1/sample_overview"),
    )
    args = parser.parse_args()
    paths = DataPaths.from_environment()
    overview = args.overview or (
        paths.data_root / "canonical" / "phenotypes" / "sample_overview_hub_geo_v1.parquet"
    )
    frame = pd.read_parquet(overview)
    # Re-apply widened HM450 band for any leftover 400-550k probe rows.
    if "n_probes" in frame.columns and "platform_id" in frame.columns:
        need = frame["platform_id"].isna() & frame["n_probes"].notna()
        for idx in frame.index[need]:
            plat = platform_from_n_probes(int(frame.at[idx, "n_probes"]))
            if plat in CATALOG_PLATFORMS:
                frame.at[idx, "platform_id"] = plat
                frame.at[idx, "platform_source"] = "ewas_db_n_probes"
        frame.to_parquet(overview, index=False)
        csv_path = overview.parent / "sample_overview_hub_geo_v1.csv.gz"
        frame.to_csv(csv_path, index=False, compression="gzip")

    summary = build_summary(frame, data_root=paths.data_root)
    report_dir = args.report_dir
    report_dir.mkdir(parents=True, exist_ok=True)
    fig_paths = write_figures(summary, report_dir / "figures")
    summary["figures"] = fig_paths
    pdf_path = report_dir / "sample_overview_dashboard.pdf"
    write_pdf_dashboard(summary, pdf_path)
    summary["pdf"] = str(pdf_path)
    write_json(report_dir / "summary.json", summary)
    write_markdown(summary, report_dir / "summary.md")
    print(  # noqa: T201
        json.dumps(
            {
                "report_dir": str(report_dir),
                "n_samples": summary["n_samples"],
                "pdf": str(pdf_path),
                "platform_top": summary["platform_top"],
                "species_status": summary["species_status"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
