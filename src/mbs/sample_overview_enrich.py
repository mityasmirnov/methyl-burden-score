"""Platform / species enrichment helpers for the sample overview census."""

from __future__ import annotations

import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

from mbs.geo_metadata import (
    SPECIES_HUMAN,
    SPECIES_NON_HUMAN,
    SPECIES_UNKNOWN,
    catalog_platform_from_gpl,
    classify_species,
)
from mbs.platform_id import normalize_platform

# Human consortium / project namespaces (no organism column in catalog).
HUMAN_PROJECT_PREFIXES: tuple[str, ...] = (
    "TCGA",
    "CPTAC",
    "ENCODE",
    "TARGET",
    "CGCI",
    "HTMCP",
    "HipSci",
    "HPSI",
)

CATALOG_PLATFORMS = frozenset({"HM450", "EPIC", "EPICv2"})


def count_text_lines(path: Path) -> int:
    """Count newline-terminated rows in a probe×beta text file (binary scan)."""
    n = 0
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            n += chunk.count(b"\n")
    return n


def platform_from_n_probes(n_probes: int) -> str | None:
    """Map Illumina-like probe counts to catalog platform_id."""
    if n_probes <= 0:
        return None
    if 400_000 <= n_probes <= 550_000:
        return "HM450"
    if 800_000 <= n_probes <= 900_000:
        return "EPIC"
    if 910_000 <= n_probes <= 980_000:
        return "EPICv2"
    return None


def platform_from_sample_id_hint(sample_id: str) -> str | None:
    """EWAS_db EPICv2 files often use a ``_935k`` sample_id suffix."""
    text = sample_id.strip().lower()
    if text.endswith("_935k") or "_935k" in text:
        return "EPICv2"
    return None


def datahub_platform_from_metadata_json(raw: object) -> str | None:
    if raw is None or (isinstance(raw, float) and pd.isna(raw)):
        return None
    try:
        payload = json.loads(raw) if isinstance(raw, str) else raw
    except (TypeError, json.JSONDecodeError):
        return None
    if not isinstance(payload, dict):
        return None
    datahub = payload.get("datahub")
    if not isinstance(datahub, dict):
        return None
    return normalize_platform(datahub.get("catalog_platform_id") or datahub.get("platform"))


def species_from_series_taxon(taxon: object) -> tuple[str | None, str | None]:
    """Return (species_status, organism) only for single-species series taxa.

    Mixed ``A; B`` series stay unresolved at sample level (need SOFT).
    """
    if taxon is None or (isinstance(taxon, float) and pd.isna(taxon)):
        return None, None
    text = str(taxon).strip()
    if not text:
        return None, None
    parts = [p.strip() for p in text.split(";") if p.strip()]
    if len(parts) != 1:
        return None, None
    organism, _tax, status = classify_species(parts[0], None)
    if status == SPECIES_UNKNOWN:
        return None, None
    return status, organism


def human_project_namespace(study_id: object) -> bool:
    if study_id is None or (isinstance(study_id, float) and pd.isna(study_id)):
        return False
    text = str(study_id).strip()
    if not text:
        return False
    for prefix in HUMAN_PROJECT_PREFIXES:
        if text == prefix or text.startswith((prefix, f"{prefix}-")):
            return True
    return False


def ewas_db_sample_path(ewas_root: Path, study_id: str, sample_id: str) -> Path | None:
    """Resolve ``EWAS_db/{study}/{sample}.txt`` when present."""
    path = ewas_root / str(study_id) / f"{sample_id}.txt"
    if path.is_file():
        return path
    # CPTAC / UUID stems sometimes omit suffixes already in sample_id.
    base = sample_id
    if base.endswith("_935k"):
        alt = ewas_root / str(study_id) / f"{base}.txt"
        if alt.is_file():
            return alt
    return path if path.is_file() else None


def infer_platforms_for_rows(
    frame: pd.DataFrame,
    *,
    data_root: Path,
    max_workers: int = 8,
) -> pd.DataFrame:
    """Fill ``platform_id`` / provenance / ``n_probes`` where missing.

    Priority: existing → DataHub metadata_json → GEO GPL map → ``_935k`` hint →
    EWAS_db probe-count bands.
    """
    out = frame.copy()
    if "platform_source" not in out.columns:
        out["platform_source"] = pd.Series([None] * len(out), dtype=object)
    if "n_probes" not in out.columns:
        out["n_probes"] = pd.Series([pd.NA] * len(out), dtype="Int64")

    # Mark already-known platforms.
    known = out["platform_id"].notna() & (
        out["platform_id"].astype(str).str.strip().isin(CATALOG_PLATFORMS)
        | out["platform_id"].astype(str).str.strip().isin({"HM450", "EPIC", "EPICv2"})
    )
    out.loc[known & out["platform_source"].isna(), "platform_source"] = "explicit"

    # DataHub sample-level platform from catalog metadata_json.
    if "metadata_json" in out.columns:
        need = ~known
        if bool(need.any()):
            dh = out.loc[need, "metadata_json"].map(datahub_platform_from_metadata_json)
            hit = dh.isin(CATALOG_PLATFORMS)
            idxs = dh.index[hit]
            out.loc[idxs, "platform_id"] = dh.loc[idxs].to_numpy()
            out.loc[idxs, "platform_source"] = "datahub_metadata_json"
        known = out["platform_id"].notna()

    # GEO GPL → catalog.
    if "geo_gpl_id" in out.columns:
        need = ~known
        if bool(need.any()):
            gpl = out.loc[need, "geo_gpl_id"].map(
                lambda x: catalog_platform_from_gpl(None if pd.isna(x) else str(x))
            )
            hit = gpl.isin(CATALOG_PLATFORMS)
            idxs = gpl.index[hit]
            out.loc[idxs, "platform_id"] = gpl.loc[idxs].to_numpy()
            out.loc[idxs, "platform_source"] = "geo_gpl"
        known = out["platform_id"].notna()

    # Filename / sample_id EPICv2 hint.
    need = ~known
    if bool(need.any()):
        hint = out.loc[need, "sample_id"].astype(str).map(platform_from_sample_id_hint)
        hit = hint.isin(CATALOG_PLATFORMS)
        idxs = hint.index[hit]
        out.loc[idxs, "platform_id"] = hint.loc[idxs].to_numpy()
        out.loc[idxs, "platform_source"] = "sample_id_935k_hint"
    known = out["platform_id"].notna()

    # Probe-count inference from EWAS_db beta files.
    ewas_root = data_root / "raw" / "ewas_datahub" / "EWAS_db"
    need_idx = list(out.index[~known])
    jobs: list[tuple[object, Path]] = []
    for idx in need_idx:
        study_id = str(out.at[idx, "study_id"])
        sample_id = str(out.at[idx, "sample_id"])
        path = ewas_db_sample_path(ewas_root, study_id, sample_id)
        if path is not None and path.is_file():
            jobs.append((idx, path))

    def _job(item: tuple[object, Path]) -> tuple[object, int | None, str | None]:
        idx, path = item
        try:
            n = count_text_lines(path)
        except OSError:
            return idx, None, None
        return idx, n, platform_from_n_probes(n)

    if jobs:
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            futures = [pool.submit(_job, job) for job in jobs]
            for fut in as_completed(futures):
                idx, n_probes, plat = fut.result()
                if n_probes is not None:
                    out.at[idx, "n_probes"] = int(n_probes)
                if plat in CATALOG_PLATFORMS:
                    out.at[idx, "platform_id"] = plat
                    out.at[idx, "platform_source"] = "ewas_db_n_probes"

    still = out["platform_id"].isna()
    out.loc[still, "platform_source"] = out.loc[still, "platform_source"].fillna("unresolved")
    return out


def fill_species_for_rows(frame: pd.DataFrame) -> pd.DataFrame:
    """Ensure every row has ``species_status`` + ``species_source`` provenance."""
    out = frame.copy()
    if "species_source" not in out.columns:
        out["species_source"] = pd.Series([None] * len(out), dtype=object)

    # Normalize empty / literal "missing" from prior overview builds.
    if "species_status" in out.columns:
        status = out["species_status"]
        blank = status.isna() | (status.astype(str).str.strip() == "") | (
            status.astype(str).str.strip().str.lower() == "missing"
        )
        out.loc[blank, "species_status"] = None
    else:
        out["species_status"] = None

    # 1) GEO sample already classified.
    geo_known = out["species_status"].isin([SPECIES_HUMAN, SPECIES_NON_HUMAN])
    out.loc[geo_known & out["species_source"].isna(), "species_source"] = "geo_soft"

    # 2) Single-species GEO series taxon (vectorized apply on unresolved only).
    if "series_taxon" in out.columns:
        need = out["species_status"].isna()
        if bool(need.any()):
            parsed = out.loc[need, "series_taxon"].map(species_from_series_taxon)
            series_status = parsed.map(lambda x: x[0])
            series_org = parsed.map(lambda x: x[1])
            hit = series_status.notna()
            idxs = series_status.index[hit]
            out.loc[idxs, "species_status"] = series_status.loc[idxs].to_numpy()
            out.loc[idxs, "species_source"] = "geo_series_taxon"
            if "organism" in out.columns:
                org_blank = out.loc[idxs, "organism"].isna() | (
                    out.loc[idxs, "organism"].astype(str).str.strip() == ""
                )
                fill_org = idxs[org_blank.to_numpy()]
                out.loc[fill_org, "organism"] = series_org.loc[fill_org].to_numpy()
            human_idxs = idxs[series_status.loc[idxs].to_numpy() == SPECIES_HUMAN]
            if "taxon_id" in out.columns:
                out.loc[human_idxs, "taxon_id"] = out.loc[human_idxs, "taxon_id"].fillna(9606)
            else:
                out.loc[human_idxs, "taxon_id"] = 9606

    # 3) Hub baseline → human (Hub∩GEO has zero non-human empirically).
    need = out["species_status"].isna()
    hub = need & out["in_hub_baseline"].fillna(False).astype(bool)
    out.loc[hub, "species_status"] = SPECIES_HUMAN
    out.loc[hub, "species_source"] = "hub_baseline_assumed_human"
    if "taxon_id" in out.columns:
        out.loc[hub, "taxon_id"] = out.loc[hub, "taxon_id"].fillna(9606)
    else:
        out.loc[hub, "taxon_id"] = 9606

    # 4) Human consortium namespaces.
    need = out["species_status"].isna()
    if bool(need.any()):
        ns_hit = out.loc[need, "study_id"].map(human_project_namespace)
        ns_idxs = ns_hit.index[ns_hit.to_numpy()]
        out.loc[ns_idxs, "species_status"] = SPECIES_HUMAN
        out.loc[ns_idxs, "species_source"] = "project_namespace_human"
        if "taxon_id" not in out.columns:
            out["taxon_id"] = pd.Series([pd.NA] * len(out), dtype="Int64")
        out.loc[ns_idxs, "taxon_id"] = out.loc[ns_idxs, "taxon_id"].fillna(9606)

    # 5) Remaining EWAS_db catalog rows: assume human until SOFT overrides.
    need = out["species_status"].isna()
    out.loc[need, "species_status"] = SPECIES_HUMAN
    out.loc[need, "species_source"] = "ewas_datahub_assumed_human"
    if "taxon_id" in out.columns:
        out.loc[need, "taxon_id"] = out.loc[need, "taxon_id"].fillna(9606)
    else:
        out.loc[need, "taxon_id"] = 9606

    if out["species_status"].isna().any():
        raise RuntimeError("species_status still null after fill rules")
    return out
