#!/usr/bin/env python3
"""Generate PNG figures for pipeline docs / inspection reports.

Requires: ``uv sync --extra analysis`` (matplotlib).
Reads ``annotation_coverage_v1``, ``annotation_graph_v1``, ``raw_inventory``,
and ``ewas_metadata_structure`` JSON; writes under each report's ``figures/``.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROLE_ORDER = (
    "promoter_core",
    "promoter_proximal",
    "five_prime",
    "gene_body",
    "three_prime",
)
PLATFORM_ORDER = ("HM450", "EPIC", "EPICv2")


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _save(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def fig_platform_assigned(cov: dict, out: Path) -> None:
    plats = {p["platform_id"]: p for p in cov["probe_level"]["platforms"]}
    labels = [p for p in PLATFORM_ORDER if p in plats]
    assigned = [plats[p]["pct_mapped_assigned"] for p in labels]
    unassigned = [plats[p]["pct_mapped_unassigned"] for p in labels]
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.bar(x, assigned, label="Mapped → typed region", color="#2a6f97")
    ax.bar(x, unassigned, bottom=assigned, label="Mapped unassigned", color="#ee9b00")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("% of mapped probes")
    ax.set_ylim(0, 100)
    ax.set_title("Per-array region assignment (mapped probes)")
    ax.axhline(70.34, color="#6c757d", ls="--", lw=1)
    ax.legend(loc="lower right", frameon=False)
    _save(fig, out)


def fig_loci_roles(cov: dict, out: Path) -> None:
    loc = cov["locus_level"]
    counts = [loc["loci_by_role"][r] for r in ROLE_ORDER]
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    ax.barh(list(ROLE_ORDER)[::-1], counts[::-1], color="#264653")
    ax.set_xlabel("Unique loci with ≥1 edge of that role")
    ax.set_title("Canonical loci by regulatory role")
    _save(fig, out)


def fig_locus_pie(cov: dict, out: Path) -> None:
    loc = cov["locus_level"]
    fig, ax = plt.subplots(figsize=(5.5, 5.0))
    ax.pie(
        [loc["n_assigned_loci"], loc["n_unassigned_loci"]],
        labels=[
            f"Assigned\n{loc['n_assigned_loci']:,}\n({loc['pct_assigned']:.1f}%)",
            f"Unassigned\n{loc['n_unassigned_loci']:,}\n({loc['pct_unassigned']:.1f}%)",
        ],
        colors=["#2a9d8f", "#e9c46a"],
        startangle=90,
        wedgeprops={"linewidth": 1, "edgecolor": "white"},
    )
    ax.set_title("Canonical loci: typed region vs unassigned")
    _save(fig, out)


def fig_cgi_context(graph: dict, out: Path) -> None:
    cgi = graph.get("cpg_context_counts") or {}
    order = [
        "island",
        "north_shore",
        "south_shore",
        "north_shelf",
        "south_shelf",
        "open_sea",
    ]
    labels = [k for k in order if k in cgi]
    vals = [cgi[k] for k in labels]
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    ax.bar(labels, vals, color="#457b9d")
    ax.set_ylabel("Loci")
    ax.set_title("CpG-island context (orthogonal locus feature)")
    ax.tick_params(axis="x", rotation=25)
    _save(fig, out)


def fig_hub_pack_bytes(inv: dict, out: Path) -> None:
    packs = inv["hub_profile_packs"]
    names = [p["family"] for p in packs]
    adv = [float(p["advertised_gb"]) for p in packs]
    ondisk = [(p["bytes"] / (1024**3)) if p.get("bytes") is not None else 0.0 for p in packs]
    x = np.arange(len(names))
    w = 0.38
    fig, ax = plt.subplots(figsize=(9.0, 4.4))
    ax.bar(x - w / 2, adv, w, label="Advertised GB", color="#8d99ae")
    ax.bar(x + w / 2, ondisk, w, label="On-disk GiB", color="#d62828")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=30, ha="right")
    ax.set_ylabel("Size")
    ax.set_title("Hub profile packs: advertised vs on-disk")
    ax.legend(frameon=False)
    # Mark disease bad zip
    for i, p in enumerate(packs):
        if p.get("status") != "ok":
            ax.annotate(
                p.get("status", "?"),
                (i + w / 2, ondisk[i]),
                textcoords="offset points",
                xytext=(0, 6),
                ha="center",
                fontsize=8,
                color="#d62828",
            )
    _save(fig, out)


def fig_hub_samples(meta: dict, out: Path) -> None:
    packs = meta.get("sample_packs") or []
    rows = sorted(
        [
            (p["family"], int(p.get("n_sample_id") or p.get("n_rows") or 0))
            for p in packs
            if p.get("exists", True)
        ],
        key=lambda t: -t[1],
    )
    if not rows:
        return
    names, counts = zip(*rows, strict=True)
    fig, ax = plt.subplots(figsize=(8.5, 4.4))
    ax.barh(list(names)[::-1], list(counts)[::-1], color="#1d3557")
    ax.set_xlabel("Unique sample_id (GSM)")
    ax.set_title("Hub sample-info unique samples by family")
    _save(fig, out)


def fig_raw_trees(inv: dict, out: Path) -> None:
    trees = inv.get("trees") or {}
    keys = ["ewas_datahub", "cpgcorpus", "manifests", "ewas_atlas"]
    labels = []
    vals = []
    for k in keys:
        t = trees.get(k) or {}
        b = t.get("bytes")
        if b is None:
            continue
        labels.append(k)
        vals.append(b / (1024**3))
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    ax.bar(labels, vals, color="#023e8a")
    ax.set_ylabel("GiB")
    ax.set_title("On-disk size by raw source lane")
    ax.tick_params(axis="x", rotation=15)
    _save(fig, out)


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    insp = root / "reports" / "inspection"
    cov = _load(insp / "annotation_coverage_v1" / "summary.json")
    graph = _load(insp / "annotation_graph_v1" / "summary.json")
    inv = _load(insp / "raw_inventory" / "summary.json")
    meta = _load(insp / "ewas_metadata_structure" / "summary.json")

    cov_fig = insp / "annotation_coverage_v1" / "figures"
    inv_fig = insp / "raw_inventory" / "figures"
    fig_platform_assigned(cov, cov_fig / "platform_assigned_vs_unassigned.png")
    fig_loci_roles(cov, cov_fig / "loci_by_role.png")
    fig_locus_pie(cov, cov_fig / "locus_assigned_pie.png")
    fig_cgi_context(graph, cov_fig / "cpg_island_context.png")
    fig_hub_pack_bytes(inv, inv_fig / "hub_pack_sizes.png")
    fig_hub_samples(meta, inv_fig / "hub_sample_counts.png")
    fig_raw_trees(inv, inv_fig / "raw_tree_sizes.png")
    print(  # noqa: T201
        json.dumps(
            {
                "annotation_coverage_figures": str(cov_fig),
                "raw_inventory_figures": str(inv_fig),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
