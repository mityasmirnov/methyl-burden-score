"""Map Hub probe IDs to canonical GRCh38 locus identifiers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass(frozen=True, slots=True)
class ProbeLocusMap:
    """Ordered study locus columns derived from observed Hub probes."""

    locus_ids: np.ndarray  # uint64
    canonical_keys: np.ndarray  # object str
    probe_ids: np.ndarray  # object str — probe used for this column
    unmapped_probe_ids: tuple[str, ...]
    n_observed_probes: int
    n_mapped_probes: int
    n_collapsed_probes: int
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
) -> ProbeLocusMap:
    """Select one primary probe per locus among observed Hub probes.

    Unmapped probes are reported and excluded from matrix columns. Samples are
    never dropped here.
    """
    observed = pd.Series(pd.unique(pd.Series(observed_probe_ids).astype(str)), dtype="string")
    n_observed = len(observed)
    hit = edges.loc[edges["probe_id"].isin(observed.tolist())].copy()
    mapped_probes = set(hit["probe_id"].astype(str))
    unmapped = tuple(sorted(p for p in observed if p not in mapped_probes))

    if hit.empty:
        raise ValueError(
            f"none of {n_observed} observed probes map to loci for platform {platform_id}"
        )

    # One column per locus: prefer lexicographically first probe_id when several
    # observed probes map to the same locus.
    hit = hit.sort_values(["locus_id", "probe_id"], kind="mergesort")
    before = len(hit)
    hit = hit.drop_duplicates(subset=["locus_id"], keep="first").reset_index(drop=True)
    n_collapsed = before - len(hit)
    hit = hit.sort_values(["locus_id"], kind="mergesort").reset_index(drop=True)

    return ProbeLocusMap(
        locus_ids=hit["locus_id"].to_numpy(dtype=np.uint64, copy=True),
        canonical_keys=hit["canonical_key"].astype(str).to_numpy(dtype=object, copy=True),
        probe_ids=hit["probe_id"].astype(str).to_numpy(dtype=object, copy=True),
        unmapped_probe_ids=unmapped,
        n_observed_probes=n_observed,
        n_mapped_probes=len(mapped_probes),
        n_collapsed_probes=int(n_collapsed),
        platform_id=platform_id,
    )
