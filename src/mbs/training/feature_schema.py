"""Explicit CpG feature column layouts for orientation anchor and assembly."""

from __future__ import annotations

FLAT_STANDARD = "flat_standard"
FLAT_REGION = "flat_region"

# flat_standard: beta, [M], [z], static..., static_present, [norm_present]
# flat_region:   M, gene_roles, cpg_context, regulatory, 3 presence flags, observed


def m_column_index(*, feature_schema: str, include_m_value: bool) -> int | None:
    """Return the M-value column index for ragged edge features, or None."""
    if not include_m_value:
        return None
    if feature_schema == FLAT_REGION:
        return 0
    if feature_schema == FLAT_STANDARD:
        return 1
    raise ValueError(f"unknown feature_schema: {feature_schema!r}")


def observed_column_index(*, feature_schema: str) -> int | None:
    """Return observed-flag column index when present."""
    if feature_schema == FLAT_REGION:
        return -1
    return None
