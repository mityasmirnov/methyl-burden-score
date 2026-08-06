"""Canonical annotation graph builders (Stage 0).

Reads InfiniumAnnotation and GENCODE as files only — never imports vendor packages.
"""

from __future__ import annotations

from mbs.annotation.build import build_annotation_graph
from mbs.annotation.probe_ids import core_probe_id, normalize_chromosome

__all__ = [
    "build_annotation_graph",
    "core_probe_id",
    "normalize_chromosome",
]
