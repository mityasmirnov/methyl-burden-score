#!/usr/bin/env python3
"""Audit fold seed-panel manifests for provenance / sparsity / overlap.

Writes ``reports/inspection/stage0_7g_prime_seed_mask/panel_audit.md`` (+ JSON).
Fails loudly if ``graph_content_hash`` is missing or every trait selected the
full prefilter width (not demonstrably sparse).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PANEL = (
    ROOT
    / "reports/inspection/stage0_7g_prime_seed_mask/seed_panels/fold_0/seed_panel.json"
)
OUT_DIR = ROOT / "reports/inspection/stage0_7g_prime_seed_mask"


def _audit(panel: dict[str, Any]) -> dict[str, Any]:
    traits = panel.get("traits") or {}
    issues: list[str] = []
    gch = panel.get("graph_content_hash")
    if not gch:
        issues.append("graph_content_hash is null/missing (provenance blocker)")
    trait_rows: list[dict[str, Any]] = []
    for name, t in sorted(traits.items()):
        n_pre = int(t.get("n_cols_prefiltered") or 0)
        n_discovery = int(
            t.get("n_discovery_cpgs")
            or t.get("n_seed_cpgs_after_stability")
            or t.get("n_seed_cpgs")
            or t.get("n_selected_cols")
            or 0
        )
        n_pass = int(t.get("n_passing_min_frequency") or n_discovery)
        sparse = bool(t.get("sparsity_ok", n_discovery < n_pre and n_discovery > 0))
        if n_pre > 0 and n_discovery >= n_pre:
            issues.append(
                f"{name}: discovery selected {n_discovery} CpGs == prefilter width "
                f"{n_pre} — not demonstrably sparse"
            )
        if t.get("strength_cap") is not None and float(t.get("strength_cap") or 0) > 1e3:
            issues.append(
                f"{name}: strength_cap={t.get('strength_cap')} looks unscaled / nonconverged"
            )
        trait_rows.append(
            {
                "trait": name,
                "n_cols_input": t.get("n_cols_input"),
                "n_cols_prefiltered": n_pre,
                "n_zero_variance_dropped": t.get("n_zero_variance_dropped"),
                "n_discovery_cpgs": n_discovery,
                "n_seed_cpgs_after_stability": n_discovery,
                "n_seed_genes": t.get("n_seed_genes", t.get("n_genes_actual")),
                "n_expanded_gene_cpg_edges": t.get(
                    "n_expanded_gene_cpg_edges", t.get("n_enriched_locus_rows")
                ),
                "n_unique_expanded_gene_cpgs": t.get("n_unique_expanded_gene_cpgs"),
                "n_multigene_cpgs": t.get("n_multigene_cpgs"),
                "seed_fraction_of_expanded": t.get("seed_fraction_of_expanded"),
                "n_passing_min_frequency": n_pass,
                "fallback_to_top_freq": t.get("fallback_to_top_freq"),
                "max_stability_seeds": t.get("max_stability_seeds"),
                "sparsity_ok": sparse,
                "strength_cap": t.get("strength_cap"),
                "strength_cap_warning": t.get("strength_cap_warning"),
                "coefs_numerically_usable": t.get("coefs_numerically_usable"),
                "standardization": t.get("standardization"),
                "n_fits_attempted": t.get("n_fits_attempted"),
                "n_fits_converged": t.get("n_fits_converged"),
                "n_fits_nonconverged": t.get("n_fits_nonconverged"),
                "frequency_quantiles": t.get("frequency_quantiles"),
                "n_autosome_only_genes": t.get("n_autosome_only_genes"),
                "n_sex_chrom_seed_cpgs": t.get("n_sex_chrom_seed_cpgs"),
            }
        )
    if "sex_autosome" not in traits:
        issues.append("sex_autosome control trait missing")
    elif int((traits.get("sex_autosome") or {}).get("n_sex_chrom_seed_cpgs") or 0) != 0:
        issues.append("sex_autosome admits sex-chrom seed CpGs")
    return {
        "graph_content_hash": gch,
        "panel_hash": panel.get("panel_hash"),
        "fold_id": panel.get("fold_id"),
        "overlap": panel.get("overlap"),
        "traits": trait_rows,
        "issues": issues,
        "g3_matched_random_quality": panel.get("g3_matched_random_quality"),
        "ok_for_seed_mask_gpu": len(issues) == 0,
    }


def _md(audit: dict[str, Any]) -> str:
    lines = [
        "# Seed-panel audit (fold-0)",
        "",
        f"- panel_hash: `{audit.get('panel_hash')}`",
        f"- graph_content_hash: `{audit.get('graph_content_hash')}`",
        f"- ok_for_seed_mask_gpu: **{audit.get('ok_for_seed_mask_gpu')}**",
        "",
        "## Issues",
        "",
    ]
    issues = audit.get("issues") or []
    if not issues:
        lines.append("- none")
    else:
        lines.extend(f"- {x}" for x in issues)
    lines.extend(
        [
            "",
            "## Per-trait selection (ADR 0012: discovery != G2 input)",
            "",
            (
                "| trait | prefilter | discovery | genes | expanded | "
                "seed frac | sparsity_ok | strength_cap |"
            ),
            "|---|---:|---:|---:|---:|---:|:---:|---:|",
        ]
    )
    for t in audit.get("traits") or []:
        lines.append(
            f"| {t['trait']} | {t['n_cols_prefiltered']} | "
            f"{t.get('n_discovery_cpgs', t.get('n_seed_cpgs_after_stability'))} | "
            f"{t.get('n_seed_genes')} | {t.get('n_unique_expanded_gene_cpgs')} | "
            f"{t.get('seed_fraction_of_expanded')} | {t['sparsity_ok']} | "
            f"{t.get('strength_cap')} |"
        )
    ov = audit.get("overlap") or {}
    if ov:
        lines.extend(
            [
                "",
                "## Overlap",
                "",
                f"- traits: `{ov.get('traits')}`",
                f"- gene union: {ov.get('gene_union_size')}",
                f"- gene pairwise: `{ov.get('gene_pairwise_overlap')}`",
                f"- CpG union (expanded): {ov.get('cpg_union_size')}",
                f"- CpG pairwise: `{ov.get('cpg_pairwise_overlap')}`",
                f"- seed fraction of expanded: `{ov.get('seed_fraction_of_expanded')}`",
                f"- genes with only one seed CpG: `{ov.get('genes_with_only_one_seed_cpg')}`",
                f"- multi-gene CpG count: `{ov.get('multi_gene_cpg_count')}`",
            ]
        )
    g3 = audit.get("g3_matched_random_quality")
    if g3:
        lines.extend(["", "## G3 matched-random quality", "", f"`{g3}`"])
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=Path, default=DEFAULT_PANEL)
    args = parser.parse_args()
    panel = json.loads(args.panel.read_text(encoding="utf-8"))
    audit = _audit(panel)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "panel_audit.json").write_text(
        json.dumps(audit, indent=2, default=str) + "\n", encoding="utf-8"
    )
    (OUT_DIR / "panel_audit.md").write_text(_md(audit), encoding="utf-8")
    print(f"wrote {OUT_DIR / 'panel_audit.md'}")
    print(f"ok_for_seed_mask_gpu={audit['ok_for_seed_mask_gpu']} issues={len(audit['issues'])}")
    for issue in audit["issues"]:
        print(f"  - {issue}")


if __name__ == "__main__":
    main()
