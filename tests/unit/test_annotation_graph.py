"""Unit tests for Stage 0 annotation graph builders."""

from __future__ import annotations

import json
from pathlib import Path
from uuid import uuid4

import pandas as pd
import pytest

from mbs.annotation.build import build_annotation_graph
from mbs.annotation.cgi_context import annotate_cpg_context, load_cpg_islands
from mbs.annotation.gencode_regions import build_gencode_regions
from mbs.annotation.locus_registry import build_locus_registry
from mbs.annotation.manifest import validate_graph_manifest
from mbs.annotation.map_loci import map_loci_to_regions, write_regions_bed
from mbs.annotation.probe_ids import core_probe_id, normalize_chromosome


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


@pytest.fixture
def workspace(tmp_path_factory: pytest.TempPathFactory, monkeypatch: pytest.MonkeyPatch) -> Path:
    repo = _repo_root()
    base = repo / "scratch" / "pytest"
    base.mkdir(parents=True, exist_ok=True)
    # Prefer pytest tmp under /data project scratch when available.
    root = base / f"ann-{uuid4().hex}"
    root.mkdir()
    monkeypatch.setenv("MBS_ROOT", str(repo))
    monkeypatch.setenv("MBS_DATA_ROOT", str(root / "data"))
    return root


def _write_gtf(path: Path) -> None:
    # Gene A on + strand: TSS=1000; Gene B on - strand overlapping body of A.
    lines = [
        'chr1\tHAVANA\tgene\t1000\t5000\t.\t+\t.\tgene_id "ENSG00000000001.1"; gene_name "GENEA"; gene_type "protein_coding";',
        'chr1\tHAVANA\ttranscript\t1000\t5000\t.\t+\t.\tgene_id "ENSG00000000001.1"; transcript_id "ENST1.1"; gene_type "protein_coding";',
        'chr1\tHAVANA\texon\t1000\t1500\t.\t+\t.\tgene_id "ENSG00000000001.1"; transcript_id "ENST1.1"; exon_number "1"; gene_type "protein_coding";',
        'chr1\tHAVANA\tfive_prime_utr\t1000\t1100\t.\t+\t.\tgene_id "ENSG00000000001.1"; transcript_id "ENST1.1"; gene_type "protein_coding";',
        'chr1\tHAVANA\texon\t2000\t4000\t.\t+\t.\tgene_id "ENSG00000000001.1"; transcript_id "ENST1.1"; exon_number "2"; gene_type "protein_coding";',
        'chr1\tHAVANA\tthree_prime_utr\t4500\t5000\t.\t+\t.\tgene_id "ENSG00000000001.1"; transcript_id "ENST1.1"; gene_type "protein_coding";',
        'chr1\tHAVANA\tgene\t3000\t4200\t.\t-\t.\tgene_id "ENSG00000000002.1"; gene_name "GENEB"; gene_type "protein_coding";',
        'chr1\tHAVANA\ttranscript\t3000\t4200\t.\t-\t.\tgene_id "ENSG00000000002.1"; transcript_id "ENST2.1"; gene_type "protein_coding";',
        'chr1\tHAVANA\texon\t4000\t4200\t.\t-\t.\tgene_id "ENSG00000000002.1"; transcript_id "ENST2.1"; exon_number "1"; gene_type "protein_coding";',
        # Non-coding gene should be ignored.
        'chr1\tHAVANA\tgene\t10\t20\t.\t+\t.\tgene_id "ENSG00000000099.1"; gene_name "LINC"; gene_type "lncRNA";',
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _synthetic_probes() -> pd.DataFrame:
    # Positions chosen to hit each role + intergenic + multi-gene overlap.
    rows = [
        # promoter_core of GENEA (TSS 1000 ± 200)
        ("cgA_core", "HM450", "chr1", 1050),
        # promoter_proximal of GENEA (TSS-1500 .. TSS-200) →  -500..800 → pos 500
        ("cgA_prox", "HM450", "chr1", 500),
        # five_prime / first exon beyond promoter_core window
        ("cgA_five", "EPIC", "chr1", 1300),
        # gene_body
        ("cgA_body", "EPIC", "chr1", 2500),
        # three_prime
        ("cgA_three", "EPICv2", "chr1", 4600),
        # overlap GENEA body and GENEB (multi-gene)
        ("cg_multi_TC21", "EPICv2", "chr1", 3500),
        # intergenic
        ("cg_inter", "HM450", "chr1", 9000),
        # unmapped
        ("cg_bad", "HM450", None, None),
    ]
    records = []
    for probe_id, platform, chrom, pos in rows:
        mapped = chrom is not None and pos is not None
        records.append(
            {
                "probe_id": probe_id,
                "platform_id": platform,
                "probe_design": "2",
                "core_probe_id": core_probe_id(probe_id),
                "chromosome": chrom,
                "position": pos,
                "mapping_status": "mapped" if mapped else "unmapped",
                "M_mapping": False,
                "M_nonuniq": False,
                "M_general": False,
                "mapQ": 60 if mapped else pd.NA,
                "strand": "+",
            }
        )
    return pd.DataFrame(records)


def test_core_probe_id_strips_epicv2_suffix() -> None:
    assert core_probe_id("cg00000029_TC21") == "cg00000029"
    assert core_probe_id("cg00000029_BC11") == "cg00000029"
    assert core_probe_id("cg00000029") == "cg00000029"


def test_normalize_chromosome() -> None:
    assert normalize_chromosome("1") == "chr1"
    assert normalize_chromosome("chr2") == "chr2"
    assert normalize_chromosome("MT") == "chrM"
    assert normalize_chromosome("NA") is None


def test_gencode_regions_and_precedence(workspace: Path) -> None:
    gtf = workspace / "tiny.gtf"
    _write_gtf(gtf)
    genes, regions = build_gencode_regions(gtf)
    assert set(genes["gene_id"]) == {"ENSG00000000001", "ENSG00000000002"}
    assert "promoter_core" in set(regions["region_type"])
    assert "promoter_proximal" in set(regions["region_type"])

    probes = _synthetic_probes()
    loci, _probes_out, edges = build_locus_registry(probes, cgi_path=None)
    assert len(loci) == 7  # mapped unique positions
    assert bool(edges["probe_id"].isin(["cg_multi_TC21"]).any())

    lr_edges, _rg = map_loci_to_regions(loci, regions)
    joined = edges.merge(lr_edges, on="locus_id", how="left").merge(
        regions[["region_id", "region_type", "gene_id"]], on="region_id", how="left"
    )

    def roles_for(probe: str) -> set[str]:
        subset = joined.loc[joined["probe_id"] == probe, "region_type"]
        return {str(v) for v in subset.dropna().tolist()}

    def genes_for(probe: str) -> set[str]:
        subset = joined.loc[joined["probe_id"] == probe, "gene_id"]
        return {str(v) for v in subset.dropna().tolist()}

    def roles_for_gene(probe: str, gene_id: str) -> set[str]:
        subset = joined.loc[
            (joined["probe_id"] == probe) & (joined["gene_id"] == gene_id),
            "region_type",
        ]
        return {str(v) for v in subset.dropna().tolist()}

    assert roles_for("cgA_core") == {"promoter_core"}
    assert roles_for("cgA_prox") == {"promoter_proximal"}
    assert roles_for("cgA_five") == {"five_prime"}
    assert roles_for("cgA_body") == {"gene_body"}
    assert roles_for_gene("cgA_three", "ENSG00000000001") == {"three_prime"}
    assert roles_for("cg_inter") == set()
    assert genes_for("cg_multi_TC21") == {"ENSG00000000001", "ENSG00000000002"}


def test_cgi_context_and_bed(workspace: Path) -> None:
    cgi = workspace / "cgi.bed"
    cgi.write_text("chr1\t1999\t2100\tisland1\n", encoding="utf-8")  # 0-based → 2000-2100 1-based
    islands = load_cpg_islands(cgi)
    assert islands.iloc[0]["start"] == 2000
    loci = pd.DataFrame(
        {
            "locus_id": pd.Series([1, 2, 3], dtype="uint64"),
            "chromosome": ["chr1", "chr1", "chr1"],
            "position": [2050, 2500, 12000],
        }
    )
    annotated = annotate_cpg_context(loci, cgi, shore_bp=2000, shelf_bp=4000)
    assert annotated.set_index("locus_id").loc[1, "cpg_context"] == "island"
    assert annotated.set_index("locus_id").loc[2, "cpg_context"] == "north_shore"
    assert annotated.set_index("locus_id").loc[3, "cpg_context"] == "open_sea"

    regions = pd.DataFrame(
        {
            "region_id": ["ENSG1:gene_body"],
            "gene_id": ["ENSG1"],
            "region_type": ["gene_body"],
            "chromosome": ["chr1"],
            "start": [100],
            "end": [200],
            "strand": ["+"],
        }
    )
    bed = workspace / "regions.bed"
    write_regions_bed(regions, bed)
    cols = bed.read_text(encoding="utf-8").strip().split("\t")
    assert cols[0] == "chr1"
    assert cols[1] == "99"  # 0-based
    assert cols[2] == "200"
    assert cols[3] == "ENSG1:gene_body"
    assert cols[6] == "ENSG1"
    assert cols[7] == "gene_body"


def test_build_annotation_graph_fixture(workspace: Path) -> None:
    repo = _repo_root()
    data_root = workspace / "data"
    gtf = workspace / "tiny.gtf"
    _write_gtf(gtf)
    genes, regions = build_gencode_regions(gtf)
    probes = _synthetic_probes()

    result = build_annotation_graph(
        project_root=repo,
        data_root=data_root,
        infinium_root=repo / "vendor" / "infinium_annotation",
        gencode_path=gtf,
        cgi_path=None,
        graph_id="graph-grch38-gencode38-five-role-v1",
        platforms=("HM450",),
        probes=probes,
        genes=genes,
        regions=regions,
    )
    graph_dir = Path(result["graph_dir"])
    ann_dir = Path(result["annotations_dir"])
    assert (ann_dir / "loci.parquet").is_file()
    assert (graph_dir / "genes.parquet").is_file()
    assert (graph_dir / "regions.bed").is_file()
    manifest = json.loads((graph_dir / "graph_manifest.json").read_text(encoding="utf-8"))
    validate_graph_manifest(manifest)
    assert manifest["graph_id"] == "graph-grch38-gencode38-five-role-v1"
    assert manifest["n_genes"] == 2
    assert result["validation_report"]["n_unassigned_loci"] >= 1
