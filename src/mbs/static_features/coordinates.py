"""Map MBS GRCh38 locus coordinates to CpGPT location keys."""

from __future__ import annotations


def mbs_locus_to_cpgpt_location(chromosome: str, position: int) -> str | None:
    """Convert a 1-based MBS cytosine locus to a CpGPT Ensembl location key.

    CpGPT human NTv2 dependency keys use 0-based positions without a ``chr``
    prefix (e.g. ``1:10847`` for MBS ``chr1:10848``). Alt / random / unlocalized
    contigs are unsupported and return ``None``.
    """
    chrom = str(chromosome).strip()
    if not chrom:
        return None
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    if not chrom or "_" in chrom:
        return None
    if chrom == "M":
        chrom = "MT"
    if int(position) < 1:
        return None
    return f"{chrom}:{int(position) - 1}"
