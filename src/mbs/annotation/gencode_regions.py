"""GENCODE → Stage 0 five-role gene regions."""

from __future__ import annotations

import gzip
import re
from collections import defaultdict
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import Any

import pandas as pd

from mbs.annotation.probe_ids import normalize_chromosome

REGION_TYPES: tuple[str, ...] = (
    "promoter_core",
    "promoter_proximal",
    "five_prime",
    "three_prime",
    "gene_body",
)
ROLE_PRECEDENCE: tuple[str, ...] = REGION_TYPES

# Quoted values: gene_id "ENSG..."; unquoted: exon_number 1;
_ATTR_RE = re.compile(r'(\w+)\s+(?:"([^"]*)"|([^\s";]+))')


def _parse_attrs(attr_field: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for key, quoted, bare in _ATTR_RE.findall(attr_field):
        out[key] = quoted if quoted else bare
    return out


def _strip_version(gene_id: str) -> str:
    return gene_id.split(".", 1)[0]


def _open_text(path: Path) -> Any:
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open("rt", encoding="utf-8")


def _iter_gtf_rows(path: Path) -> Iterator[tuple[str, str, int, int, str, dict[str, str]]]:
    with _open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _source, feature, start_s, end_s, _score, strand, _frame, attrs = parts[:9]
            yield chrom, feature, int(start_s), int(end_s), strand, _parse_attrs(attrs)


def _merge_intervals(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    ordered = sorted((int(a), int(b)) for a, b in intervals if a <= b)
    if not ordered:
        return []
    merged: list[tuple[int, int]] = [ordered[0]]
    for start, end in ordered[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end + 1:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged


def _promoter_windows(
    tss: int,
    strand: str,
    *,
    core_bp: int,
    proximal_bp: int,
) -> tuple[tuple[int, int], tuple[int, int] | None]:
    core = (max(1, tss - core_bp), tss + core_bp)
    if strand == "+":
        prox_end = tss - core_bp
        prox_start = tss - proximal_bp
        proximal = (max(1, prox_start), prox_end) if prox_end >= prox_start else None
    elif strand == "-":
        prox_start = tss + core_bp
        prox_end = tss + proximal_bp
        proximal = (prox_start, prox_end) if prox_end >= prox_start else None
    else:
        proximal = None
    return core, proximal


def _classify_utr(
    start: int,
    end: int,
    strand: str,
    cds_min: int,
    cds_max: int,
) -> str | None:
    """Map a GENCODE ``UTR`` interval to five_prime or three_prime via CDS bounds."""
    mid = (start + end) // 2
    role: str | None = None
    if strand == "+":
        if mid < cds_min:
            role = "five_prime"
        elif mid > cds_max:
            role = "three_prime"
    elif strand == "-":
        if mid > cds_max:
            role = "five_prime"
        elif mid < cds_min:
            role = "three_prime"
    return role


def build_gencode_regions(
    gtf_path: Path,
    *,
    source_version: str = "GENCODE_v38",
    gene_types: frozenset[str] | None = None,
    core_bp: int = 200,
    proximal_bp: int = 1500,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Parse GENCODE GTF into gene-level union regions for the five Stage 0 roles."""
    allowed = gene_types or frozenset({"protein_coding"})
    genes: dict[str, dict[str, object]] = {}
    role_intervals: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    cds_bounds: dict[str, list[tuple[int, int]]] = defaultdict(list)
    pending_utrs: list[tuple[str, str, int, int, str]] = []

    for chrom_raw, feature, start, end, strand, attrs in _iter_gtf_rows(gtf_path):
        chrom = normalize_chromosome(chrom_raw)
        if chrom is None:
            continue
        gene_type = attrs.get("gene_type") or attrs.get("gene_biotype") or ""
        if feature == "gene":
            if gene_type not in allowed:
                continue
            gene_id = _strip_version(attrs["gene_id"])
            genes[gene_id] = {
                "gene_id": gene_id,
                "gene_name": attrs.get("gene_name", gene_id),
                "chromosome": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                "gene_type": gene_type,
                "source_version": source_version,
            }
            role_intervals[gene_id]["gene_body"].append((start, end))
            continue

        if feature == "transcript":
            if gene_type not in allowed:
                continue
            gene_id = _strip_version(attrs["gene_id"])
            if gene_id not in genes:
                continue
            tss = start if strand == "+" else end
            core, proximal = _promoter_windows(
                tss, strand, core_bp=core_bp, proximal_bp=proximal_bp
            )
            role_intervals[gene_id]["promoter_core"].append(core)
            if proximal is not None:
                role_intervals[gene_id]["promoter_proximal"].append(proximal)
            continue

        gene_id = _strip_version(attrs.get("gene_id", ""))
        if not gene_id or gene_id not in genes:
            continue
        transcript_id = attrs.get("transcript_id", "")

        if feature == "exon":
            if attrs.get("exon_number") == "1":
                role_intervals[gene_id]["five_prime"].append((start, end))
        elif feature == "five_prime_utr":
            role_intervals[gene_id]["five_prime"].append((start, end))
        elif feature == "three_prime_utr":
            role_intervals[gene_id]["three_prime"].append((start, end))
        elif feature == "CDS" and transcript_id:
            cds_bounds[transcript_id].append((start, end))
        elif feature == "UTR" and transcript_id:
            pending_utrs.append((gene_id, transcript_id, start, end, strand))

    for gene_id, transcript_id, start, end, strand in pending_utrs:
        bounds = cds_bounds.get(transcript_id)
        if not bounds:
            continue
        cds_min = min(s for s, _e in bounds)
        cds_max = max(e for _s, e in bounds)
        role = _classify_utr(start, end, strand, cds_min, cds_max)
        if role is not None:
            role_intervals[gene_id][role].append((start, end))

    gene_rows = list(genes.values())
    region_rows: list[dict[str, object]] = []
    for gene_id, roles in role_intervals.items():
        gene = genes[gene_id]
        for region_type in REGION_TYPES:
            merged = _merge_intervals(roles.get(region_type, []))
            for index, (start, end) in enumerate(merged):
                suffix = "" if index == 0 else f"_{index + 1}"
                region_rows.append(
                    {
                        "region_id": f"{gene_id}:{region_type}{suffix}",
                        "gene_id": gene_id,
                        "region_type": region_type,
                        "chromosome": gene["chromosome"],
                        "start": start,
                        "end": end,
                        "strand": gene["strand"],
                        "source": "gencode",
                        "source_version": source_version,
                    }
                )

    genes_df = (
        pd.DataFrame(gene_rows).sort_values("gene_id", kind="mergesort").reset_index(drop=True)
    )
    regions_df = (
        pd.DataFrame(region_rows)
        .sort_values(["gene_id", "region_type", "start"], kind="mergesort")
        .reset_index(drop=True)
    )
    return genes_df, regions_df
