"""Map Hub probe IDs to canonical GRCh38 locus identifiers."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from mbs.batch import (
    ANNOTATION_STATUS_NAMES,
    ANNOTATION_STATUS_UNMAPPED,
)

RESIDUAL_CANONICAL_PREFIX = "residual:"
COLLAPSE_IDENTITY = "identity"
COLLAPSE_MEAN = "mean"
COLLAPSE_MEDIAN = "median"


def residual_canonical_key(probe_id: str) -> str:
    return f"{RESIDUAL_CANONICAL_PREFIX}{probe_id}"


def is_residual_canonical_key(canonical_key: str) -> bool:
    return str(canonical_key).startswith(RESIDUAL_CANONICAL_PREFIX)


def synthetic_residual_locus_id(probe_id: str) -> np.uint64:
    """Deterministic uint64 locus id for Illumina-unmapped residual columns."""
    digest = hashlib.sha256(f"mbs-residual:{probe_id}".encode()).digest()
    # Keep out of the low Illumina-mapped id space used by the locus registry.
    value = int.from_bytes(digest[:8], "big") | (1 << 63)
    return np.uint64(value)


def collapse_method_for_n(n_probes: int) -> str:
    """Collapse policy: identity (1), mean (2), median (≥3)."""
    if n_probes <= 1:
        return COLLAPSE_IDENTITY
    if n_probes == 2:
        return COLLAPSE_MEAN
    return COLLAPSE_MEDIAN


@dataclass(frozen=True, slots=True)
class ProbeLocusMap:
    """Ordered study locus columns derived from observed Hub probes.

    Mapped GRCh38 loci come first (one column per locus). Illumina-coordinate-
    unmapped probes are retained as residual columns afterward (never dropped).
    When several observed probes map to one locus, ``contributing_probe_ids``
    records all of them and ``collapse_method`` selects mean/median aggregation.
    ``probe_ids`` is a stable display id (lexicographic first) only.
    """

    locus_ids: np.ndarray  # uint64
    canonical_keys: np.ndarray  # object str
    probe_ids: np.ndarray  # object str — display probe for this column
    annotation_status: np.ndarray  # int8 — batch status ids (residual cols = unmapped)
    contributing_probe_ids: tuple[tuple[str, ...], ...]
    collapse_method: tuple[str, ...]
    unmapped_probe_ids: tuple[str, ...]
    n_observed_probes: int
    n_mapped_probes: int
    n_collapsed_probes: int
    n_residual_probes: int
    platform_id: str


def load_probe_locus_edges(
    annotations_dir: Path,
    *,
    platform_id: str,
) -> pd.DataFrame:
    """Load primary probe→locus edges for ``platform_id`` with canonical keys."""
    annotations_dir = annotations_dir.resolve()
    edges_path = annotations_dir / "probe_locus_edges.parquet"
    loci_path = annotations_dir / "loci.parquet"
    if not edges_path.is_file():
        raise FileNotFoundError(f"missing probe_locus_edges: {edges_path}")
    if not loci_path.is_file():
        raise FileNotFoundError(f"missing loci registry: {loci_path}")

    edges = pd.read_parquet(
        edges_path,
        columns=["probe_id", "platform_id", "locus_id", "is_primary"],
    )
    edges = edges.loc[
        (edges["platform_id"] == platform_id) & (edges["is_primary"].astype(bool))
    ].copy()
    if edges.empty:
        raise ValueError(f"no primary probe_locus_edges for platform_id={platform_id!r}")

    loci = pd.read_parquet(loci_path, columns=["locus_id", "canonical_key", "genome_build"])
    builds = set(loci["genome_build"].astype(str).unique())
    if builds != {"GRCh38"}:
        raise ValueError(f"locus registry genome_build must be GRCh38 only, got {sorted(builds)}")

    merged = edges.merge(loci, on="locus_id", how="inner", validate="many_to_one")
    if merged.empty:
        raise ValueError(f"probe edges for {platform_id} did not join any loci")
    return merged.sort_values(["probe_id", "locus_id"], kind="mergesort").reset_index(drop=True)


def build_probe_locus_map(
    observed_probe_ids: np.ndarray | list[str] | pd.Series,
    edges: pd.DataFrame,
    *,
    platform_id: str,
    retain_unmapped: bool = True,
) -> ProbeLocusMap:
    """Build one matrix column per locus among observed Hub probes.

    Illumina-coordinate-unmapped probes are retained as residual columns when
    ``retain_unmapped`` is True (Milestone 6 default). Samples are never dropped.
    Multiple observed probes at one locus are collapsed (mean/median), not dropped.
    """
    observed = pd.Series(pd.unique(pd.Series(observed_probe_ids).astype(str)), dtype="string")
    n_observed = len(observed)
    hit = edges.loc[edges["probe_id"].isin(observed.tolist())].copy()
    mapped_probes = set(hit["probe_id"].astype(str))
    unmapped = tuple(sorted(p for p in observed if p not in mapped_probes))

    if hit.empty and not (retain_unmapped and unmapped):
        raise ValueError(
            f"none of {n_observed} observed probes map to loci for platform {platform_id}"
        )

    locus_ids: list[np.uint64] = []
    canonical_keys: list[str] = []
    probe_ids: list[str] = []
    statuses: list[int] = []
    contributing: list[tuple[str, ...]] = []
    methods: list[str] = []
    n_collapsed = 0

    if not hit.empty:
        hit = hit.sort_values(["locus_id", "probe_id"], kind="mergesort")
        for locus_id, group in hit.groupby("locus_id", sort=True):
            probes = tuple(sorted(str(p) for p in group["probe_id"].tolist()))
            method = collapse_method_for_n(len(probes))
            if len(probes) > 1:
                n_collapsed += len(probes) - 1
            canonical = str(group["canonical_key"].iloc[0])
            locus_ids.append(np.uint64(locus_id))
            canonical_keys.append(canonical)
            probe_ids.append(probes[0])  # display id only
            contributing.append(probes)
            methods.append(method)
            # Regulatory status is finalized at train time; Illumina-mapped columns
            # start as mapped placeholders until graph join.
            statuses.append(0)  # ANNOTATION_STATUS_MAPPED placeholder

    n_residual = 0
    if retain_unmapped and unmapped:
        for probe_id in unmapped:
            locus_ids.append(synthetic_residual_locus_id(probe_id))
            canonical_keys.append(residual_canonical_key(probe_id))
            probe_ids.append(probe_id)
            contributing.append((probe_id,))
            methods.append(COLLAPSE_IDENTITY)
            statuses.append(ANNOTATION_STATUS_UNMAPPED)
            n_residual += 1

    if not locus_ids:
        raise ValueError(f"no matrix columns for platform {platform_id}")

    return ProbeLocusMap(
        locus_ids=np.asarray(locus_ids, dtype=np.uint64),
        canonical_keys=np.asarray(canonical_keys, dtype=object),
        probe_ids=np.asarray(probe_ids, dtype=object),
        annotation_status=np.asarray(statuses, dtype=np.int8),
        contributing_probe_ids=tuple(contributing),
        collapse_method=tuple(methods),
        unmapped_probe_ids=unmapped,
        n_observed_probes=n_observed,
        n_mapped_probes=len(mapped_probes),
        n_collapsed_probes=int(n_collapsed),
        n_residual_probes=n_residual,
        platform_id=platform_id,
    )


def annotation_status_name(status_id: int) -> str:
    if status_id < 0 or status_id >= len(ANNOTATION_STATUS_NAMES):
        raise ValueError(f"unknown annotation status id {status_id}")
    return ANNOTATION_STATUS_NAMES[status_id]
