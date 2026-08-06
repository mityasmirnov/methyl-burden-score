"""Read Zhou-lab InfiniumAnnotation tables without importing vendor code."""

from __future__ import annotations

import gzip
from collections.abc import Iterable
from pathlib import Path

import numpy as np
import pandas as pd

from mbs.annotation.probe_ids import core_probe_id, normalize_chromosome

DEFAULT_PLATFORMS: tuple[str, ...] = ("HM450", "EPIC", "EPICv2")
MASK_TAGS: tuple[str, ...] = ("M_mapping", "M_nonuniq", "M_general")


def _read_ordering(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype={"Probe_ID": "string"})


def _read_coord(path: Path) -> pd.DataFrame:
    return pd.read_csv(
        path,
        sep="\t",
        dtype={
            "CpG_chrm": "string",
            "CpG_beg": "Int64",
            "strand": "string",
            "mapQ": "Int64",
        },
    )


def read_yame_mask_bits(
    cm_path: Path,
    idx_path: Path,
    n_probes: int,
    tags: Iterable[str] = MASK_TAGS,
) -> dict[str, np.ndarray]:
    """Unpack selected YAME ``.cm`` bitsets aligned to ordering rows.

    The compressed payload is a short header followed by one packed bitset per
    mask listed in the ``.idx`` file (order preserved). Bit order is little-endian
    within each byte (numpy ``bitorder='little'``).
    """
    names = [line.split("\t")[0] for line in idx_path.read_text().splitlines() if line.strip()]
    with gzip.open(cm_path, "rb") as handle:
        raw = handle.read()
    bytes_per = (n_probes + 7) // 8
    header_size = len(raw) - len(names) * bytes_per
    if header_size < 0:
        raise ValueError(f"mask payload shorter than expected bitsets: {cm_path}")
    data = raw[header_size:]
    wanted = set(tags)
    out: dict[str, np.ndarray] = {}
    for index, name in enumerate(names):
        if name not in wanted:
            continue
        chunk = data[index * bytes_per : (index + 1) * bytes_per]
        bits = np.unpackbits(np.frombuffer(chunk, dtype=np.uint8), bitorder="little")[:n_probes]
        out[name] = bits.astype(bool, copy=False)
    missing = wanted - set(out)
    if missing:
        raise KeyError(f"mask tags not found in {idx_path}: {sorted(missing)}")
    return out


def load_platform_probes(infinium_root: Path, platform_id: str) -> pd.DataFrame:
    """Load one platform's probes with GRCh38 coords and QC mask flags."""
    plat_dir = infinium_root / platform_id
    ordering_path = plat_dir / f"{platform_id}.ordering.tsv.gz"
    coord_path = plat_dir / f"{platform_id}.hg38.coord.tsv.gz"
    mask_path = plat_dir / f"{platform_id}.hg38.mask.cm"
    idx_path = plat_dir / f"{platform_id}.hg38.mask.cm.idx"
    for path in (ordering_path, coord_path, mask_path, idx_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    ordering = _read_ordering(ordering_path)
    coord = _read_coord(coord_path)
    if len(ordering) != len(coord):
        raise ValueError(
            f"{platform_id}: ordering ({len(ordering)}) and coord ({len(coord)}) length mismatch"
        )

    masks = read_yame_mask_bits(mask_path, idx_path, n_probes=len(ordering))
    frame = pd.DataFrame(
        {
            "probe_id": ordering["Probe_ID"].astype(str),
            "platform_id": platform_id,
            "probe_design": ordering["col"].astype("string"),
            "chrom_raw": coord["CpG_chrm"],
            "cpg_beg_0based": coord["CpG_beg"],
            "strand": coord["strand"].astype("string"),
            "mapQ": coord["mapQ"],
            "M_mapping": masks["M_mapping"],
            "M_nonuniq": masks["M_nonuniq"],
            "M_general": masks["M_general"],
        }
    )
    frame["core_probe_id"] = [core_probe_id(str(pid)) for pid in frame["probe_id"]]
    frame["chromosome"] = [
        normalize_chromosome(str(value)) if pd.notna(value) else None
        for value in frame["chrom_raw"]
    ]
    # InfiniumAnnotation CpG_beg is 0-based cytosine start → store 1-based position.
    mapped = frame["chromosome"].notna() & frame["cpg_beg_0based"].notna()
    frame["position"] = pd.Series(pd.NA, index=frame.index, dtype="Int64")
    frame.loc[mapped, "position"] = frame.loc[mapped, "cpg_beg_0based"].astype("int64") + 1
    frame["mapping_status"] = np.where(mapped, "mapped", "unmapped")
    return frame


def load_infinium_probes(
    infinium_root: Path,
    platforms: Iterable[str] = DEFAULT_PLATFORMS,
) -> pd.DataFrame:
    """Concatenate platform probe tables."""
    frames = [load_platform_probes(infinium_root, platform_id) for platform_id in platforms]
    if not frames:
        raise ValueError("no platforms requested")
    return pd.concat(frames, ignore_index=True)
