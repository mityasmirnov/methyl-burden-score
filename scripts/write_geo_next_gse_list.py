#!/usr/bin/env python3
"""Pick next GEO backfill GSE list from EWAS_db on disk.

Starts from EWAS_db study dirs with ≥1 GSM ``.txt``, minus already-fetched
studies (cached batch / parquet / optional exclude file). Ranks by assay N;
when family SOFT is already cached, boosts age/tissue coverage and drops
non-human series.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from mbs.annotation.manifest import utc_now_iso, write_json
from mbs.geo_metadata import (
    build_geo_frame_from_soft,
    load_geo_frame,
    read_cached_soft,
    species_census,
)
from mbs.paths import DataPaths

_GSM_TXT = re.compile(r"^GSM[0-9]+\.txt$", re.I)
DEFAULT_EXCLUDE = Path("configs/data/geo_backfill_cached_batch_gse.txt")
DEFAULT_OUT = Path("configs/data/geo_backfill_next_gse.txt")
DEFAULT_REPORT = Path("reports/inspection/deepmat_data_v1/geo_backfill_next/summary.json")


def _load_ids(path: Path) -> set[str]:
    if not path.is_file():
        return set()
    return {
        line.strip().upper()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith("#")
    }


def _n_gsm(study_dir: Path) -> int:
    return sum(1 for f in study_dir.iterdir() if f.is_file() and _GSM_TXT.match(f.name))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=100, help="Max GSE to write (default 100)")
    parser.add_argument("--min-gsm", type=int, default=50, help="Minimum assay GSM on disk")
    parser.add_argument(
        "--exclude-file",
        type=Path,
        default=DEFAULT_EXCLUDE,
        help="Already-fetched GSE list (default: cached batch)",
    )
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    paths = DataPaths.from_environment()
    ewas_root = paths.data_root / "raw" / "ewas_datahub" / "EWAS_db"
    if not ewas_root.is_dir():
        raise SystemExit(f"missing EWAS_db root: {ewas_root}")

    exclude = _load_ids(args.exclude_file.resolve())
    # Also exclude studies already in the GEO parquet.
    geo = load_geo_frame(paths.data_root)
    if not geo.empty and "study_id" in geo.columns:
        exclude |= {
            str(s).strip().upper()
            for s in geo["study_id"].dropna().astype(str).tolist()
            if str(s).strip()
        }

    candidates: list[dict[str, object]] = []
    for study_dir in sorted(ewas_root.iterdir()):
        if not study_dir.is_dir():
            continue
        gse = study_dir.name.strip().upper()
        if not gse.startswith("GSE"):
            continue
        if gse in exclude:
            continue
        n_assay = _n_gsm(study_dir)
        if n_assay < args.min_gsm:
            continue
        soft = paths.cache_root / "geo" / gse / f"{gse}_family.soft.gz"
        row: dict[str, object] = {
            "study_id": gse,
            "n_assay_gsm": n_assay,
            "soft_cached": soft.is_file(),
            "n_age": None,
            "n_tissue_mapped": None,
            "species_human": None,
            "species_non_human": None,
            "species_unknown": None,
            "species_ok": True,
            "score": float(n_assay),
        }
        if soft.is_file():
            try:
                frame = build_geo_frame_from_soft(
                    read_cached_soft(soft),
                    fetched_at=utc_now_iso(),
                    soft_sha256="cache",
                )
                census = species_census(frame)
                row["species_human"] = census["human"]
                row["species_non_human"] = census["non_human"]
                row["species_unknown"] = census["unknown"]
                row["n_age"] = int(frame["age"].notna().sum()) if "age" in frame.columns else 0
                row["n_tissue_mapped"] = (
                    int((frame["tissue_map_status"] == "mapped").sum())
                    if "tissue_map_status" in frame.columns
                    else 0
                )
                # Drop series with any non-human GSM (quarantine whole GSE for next list).
                if int(census["non_human"]) > 0:
                    row["species_ok"] = False
                # Prefer age/tissue-rich when SOFT already known.
                row["score"] = float(n_assay) + 2.0 * float(row["n_age"] or 0) + 2.0 * float(
                    row["n_tissue_mapped"] or 0
                )
            except OSError:
                pass
        if row["species_ok"]:
            candidates.append(row)

    candidates.sort(key=lambda r: (-float(r["score"]), str(r["study_id"])))
    selected = candidates[: max(0, args.limit)]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Next GEO backfill GSE list (EWAS_db on disk, not yet fetched).",
        "# Ranked by assay GSM (+ age/tissue boost when SOFT cached).",
        f"# Generated: {utc_now_iso()}",
        f"# exclude={args.exclude_file} min_gsm={args.min_gsm} limit={args.limit}",
        "",
    ]
    lines.extend(str(r["study_id"]) for r in selected)
    lines.append("")
    args.out.write_text("\n".join(lines), encoding="utf-8")

    report = {
        "generated_at": utc_now_iso(),
        "exclude_file": str(args.exclude_file),
        "n_excluded": len(exclude),
        "min_gsm": args.min_gsm,
        "limit": args.limit,
        "n_candidates": len(candidates),
        "n_selected": len(selected),
        "out": str(args.out),
        "notes": [
            "Species gate: series with any non-human GSM (when SOFT cached) are omitted.",
            "Uncached SOFTs are ranked by assay N only; species checked at fetch time.",
            "Age/tissue richness is unknown until SOFT fetch for uncached studies.",
            "Larger crawl remains gated for training; this only picks the next fetch list.",
        ],
        "selected": selected,
        "top_held_out": candidates[args.limit : args.limit + 20],
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    write_json(args.report, report)
    md = args.report.with_suffix(".md")
    md_lines = [
        "# Next GEO backfill GSE list",
        "",
        f"- Generated: `{report['generated_at']}`",
        f"- Candidates (≥{args.min_gsm} GSM, not fetched): **{len(candidates)}**",
        f"- Selected: **{len(selected)}** → `{args.out}`",
        "",
        "| study_id | n_assay | soft | age | tissue_mapped | score |",
        "| --- | ---: | --- | ---: | ---: | ---: |",
    ]
    for r in selected[:40]:
        md_lines.append(
            f"| `{r['study_id']}` | {r['n_assay_gsm']} | {r['soft_cached']} | "
            f"{r['n_age']} | {r['n_tissue_mapped']} | {r['score']:.0f} |"
        )
    md_lines.extend(["", "## Notes", ""])
    md_lines.extend(f"- {n}" for n in report["notes"])
    md_lines.append("")
    md.write_text("\n".join(md_lines), encoding="utf-8")
    print(f"wrote {args.out} n={len(selected)}")
    print(f"wrote {args.report}")
    print(f"wrote {md}")


if __name__ == "__main__":
    main()
