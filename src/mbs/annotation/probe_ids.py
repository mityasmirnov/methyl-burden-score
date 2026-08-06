"""Probe identifier and chromosome helpers."""

from __future__ import annotations

import re

# Illumina EPICv2 technical suffixes, e.g. cg00000029_TC21 or cg00000029_BC11.
_EPICV2_SUFFIX_RE = re.compile(r"_[A-Z]{2}\d{2}$")


def core_probe_id(probe_id: str) -> str:
    """Strip EPICv2 technical suffixes; leave HM450/EPIC IDs unchanged."""
    return _EPICV2_SUFFIX_RE.sub("", probe_id.strip())


def normalize_chromosome(chrom: str) -> str | None:
    """Normalize contig labels to ``chrN`` / ``chrX`` / ``chrM`` form."""
    value = chrom.strip()
    if not value or value.upper() in {"NA", "NAN", ".", "*"}:
        return None
    if value.startswith("chr"):
        return value
    if value == "MT":
        return "chrM"
    return f"chr{value}"
