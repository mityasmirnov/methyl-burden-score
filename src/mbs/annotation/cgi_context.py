"""CpG island / shore / shelf / open-sea context for canonical loci."""

from __future__ import annotations

from pathlib import Path

import duckdb
import pandas as pd

from mbs.annotation.probe_ids import normalize_chromosome


def load_cpg_islands(cgi_path: Path) -> pd.DataFrame:
    """Load UCSC ``cpgIslandExt`` table dump or a 3+ column BED."""
    if not cgi_path.is_file():
        raise FileNotFoundError(cgi_path)

    peek = pd.read_csv(cgi_path, sep="\t", header=None, nrows=1, compression="infer")
    first = peek.iloc[0]
    # UCSC dump: bin(int), chrom, chromStart, chromEnd, ...
    # BED: chrom, start, end, ...
    if pd.api.types.is_number(first.iloc[0]) and str(first.iloc[1]).startswith("chr"):
        usecols = [1, 2, 3]
        names = ["chromosome", "start_0", "end_0"]
    else:
        usecols = [0, 1, 2]
        names = ["chromosome", "start_0", "end_0"]

    frame = pd.read_csv(
        cgi_path,
        sep="\t",
        header=None,
        compression="infer",
        usecols=usecols,
        names=names,
        dtype={"chromosome": "string", "start_0": "int64", "end_0": "int64"},
    )
    frame["chromosome"] = [
        normalize_chromosome(str(value)) if pd.notna(value) else None
        for value in frame["chromosome"]
    ]
    frame = frame.dropna(subset=["chromosome"]).copy()
    # Store 1-based inclusive interval matching locus positions.
    frame["start"] = frame["start_0"] + 1
    frame["end"] = frame["end_0"]
    return frame.loc[:, ["chromosome", "start", "end"]].reset_index(drop=True)


def annotate_cpg_context(
    loci: pd.DataFrame,
    cgi_path: Path,
    *,
    shore_bp: int = 2000,
    shelf_bp: int = 4000,
) -> pd.DataFrame:
    """Attach ``cpg_context`` using island / N-S shore / shelf / open_sea labels."""
    islands = load_cpg_islands(cgi_path)
    if loci.empty:
        out = loci.copy()
        out["cpg_context"] = pd.Series(dtype="string")
        return out

    con = duckdb.connect(database=":memory:")
    con.register("loci", loci[["locus_id", "chromosome", "position"]])
    con.register("islands", islands)
    assigned = con.execute(
        f"""
        WITH dist AS (
          SELECT
            l.locus_id,
            l.chromosome,
            l.position,
            i.start AS island_start,
            i.end AS island_end,
            CASE
              WHEN l.position BETWEEN i.start AND i.end THEN 0
              WHEN l.position < i.start THEN i.start - l.position
              ELSE l.position - i.end
            END AS abs_dist,
            CASE
              WHEN l.position BETWEEN i.start AND i.end THEN 'island'
              WHEN l.position < i.start THEN 'south'
              ELSE 'north'
            END AS side
          FROM loci l
          JOIN islands i
            ON l.chromosome = i.chromosome
           AND l.position BETWEEN (i.start - {int(shelf_bp)}) AND (i.end + {int(shelf_bp)})
        ),
        best AS (
          SELECT
            locus_id,
            side,
            abs_dist,
            ROW_NUMBER() OVER (
              PARTITION BY locus_id
              ORDER BY abs_dist ASC, side ASC
            ) AS rn
          FROM dist
        )
        SELECT
          locus_id,
          CASE
            WHEN abs_dist = 0 THEN 'island'
            WHEN abs_dist <= {int(shore_bp)} AND side = 'north' THEN 'north_shore'
            WHEN abs_dist <= {int(shore_bp)} AND side = 'south' THEN 'south_shore'
            WHEN abs_dist <= {int(shelf_bp)} AND side = 'north' THEN 'north_shelf'
            WHEN abs_dist <= {int(shelf_bp)} AND side = 'south' THEN 'south_shelf'
            ELSE 'open_sea'
          END AS cpg_context
        FROM best
        WHERE rn = 1
        """
    ).fetchdf()
    con.close()

    out = loci.drop(columns=[c for c in ("cpg_context",) if c in loci.columns]).merge(
        assigned,
        on="locus_id",
        how="left",
    )
    out["cpg_context"] = out["cpg_context"].fillna("open_sea").astype("string")
    return out
