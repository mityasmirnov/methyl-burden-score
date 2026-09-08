#!/usr/bin/env python3
"""Cheap one-hop correctness smokes (Milestone 10) — validate before scaling.

1) One-hop + seed-gene masking (G0 vs G1 on FlatDeepSetRegion). Never tested
   on one-hop (9c only covered cascade).
2) One-hop multi-seed restarts (seeds 42/43/44, fold 0).

Default: ATS fold 0, 5 epochs, mean pool, gene-linked panel, rho_hidden=64.
Reuses fold-0 seed panels from Milestone 9c when present.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from mbs.matrix.store import matrix_store_paths, open_betas_for_matrix, read_sample_index
from mbs.paths import DataPaths
from mbs.training.cascade_assign import build_cascade_assignment
from mbs.training.dev_cv import load_frozen_folds
from mbs.training.flat_region_loop import train_flat_region_on_arrays
from mbs.training.locus_gene import load_graph_tables
from mbs.training.phenotypes import load_multitask_phenotypes
from mbs.training.seed_panel import gene_mask_tensor, load_seed_panel

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "reports" / "inspection" / "stage0_7h_onehop_correctness"
LIGHT_CFG = ROOT / "configs" / "experiment" / "stage0_7g_gene_only_probe_light_mean.yaml"
SEED_PANEL_DIR = (
    ROOT / "reports" / "inspection" / "stage0_7g_prime_seed_mask" / "seed_panels" / "fold_0"
)
# N-light default after nine-pack wide screen (2026-09-08).
DEFAULT_RHO_HIDDEN = 64


def _phenotype_arrays(phenotypes: list[Any], sample_ids: list[str]) -> dict[str, np.ndarray]:
    by_id = {p.sample_id: p for p in phenotypes}
    return {
        "tissue": np.asarray([int(by_id[s].class_index) for s in sample_ids], dtype=np.int64),
        "tissue_mask": np.asarray([bool(by_id[s].tissue_mask) for s in sample_ids], dtype=bool),
        "age": np.asarray([float(by_id[s].age or 0.0) for s in sample_ids], dtype=np.float64),
        "age_mask": np.asarray([bool(by_id[s].age_mask) for s in sample_ids], dtype=bool),
        "sex": np.asarray([int(by_id[s].sex_class_index or 0) for s in sample_ids], dtype=np.int64),
        "sex_mask": np.asarray([bool(by_id[s].sex_mask) for s in sample_ids], dtype=bool),
        "study_ids": np.asarray([str(by_id[s].study_id or "NA") for s in sample_ids], dtype=object),
    }


def _split_path(paths: DataPaths, split_id: str) -> Path:
    p = paths.artifact_root / "splits" / split_id / "folds.json"
    if not p.is_file():
        raise FileNotFoundError(p)
    return p


def _load_ats_fold0(*, max_loci: int) -> dict[str, Any]:
    paths = DataPaths.from_environment()
    cfg = yaml.safe_load(LIGHT_CFG.read_text())
    pilot = cfg["pilot"]
    matrix_id = str(pilot["matrix_id"])
    graph_id = str(pilot["graph_id"])
    split_id = str(cfg.get("split_id") or pilot.get("split_id") or "hub-ats-7e-3fold-v1")
    pheno_rel = str(cfg.get("sample_phenotype_table") or pilot["sample_phenotype_table"])
    matrix_root = paths.data_root / "canonical" / "matrices" / matrix_id
    sample_index = read_sample_index(matrix_store_paths(matrix_root).sample_index_path)
    sample_ids = [str(s) for s in sample_index["sample_id"].tolist()]
    row_by_id = {sid: i for i, sid in enumerate(sample_ids)}
    phenotypes, class_names = load_multitask_phenotypes(
        paths.data_root / pheno_rel,
        sample_ids=sample_ids,
    )
    fold_pack = load_frozen_folds(_split_path(paths, split_id))
    fold0 = fold_pack["folds"][0]
    train_ids = [s for s in fold0["train_sample_ids"] if s in row_by_id]
    test_ids = [s for s in (fold0.get("external_test_sample_ids") or []) if s in row_by_id]
    if not test_ids:
        test_ids = [s for s in fold0["validation_sample_ids"] if s in row_by_id]
    train_idx = np.asarray([row_by_id[s] for s in train_ids], dtype=np.int64)
    test_idx = np.asarray([row_by_id[s] for s in test_ids], dtype=np.int64)
    arrays = _phenotype_arrays(phenotypes, sample_ids)
    loci, genes, edges = load_graph_tables(paths.data_root / "canonical" / "graphs" / graph_id)
    assignment = build_cascade_assignment(
        loci, genes, edges, gene_allocation="explicit_only"
    )
    print(
        f"[onehop-smoke] loading betas {matrix_id} cols=[:{max_loci}] …",
        flush=True,
    )
    betas = np.asarray(open_betas_for_matrix(matrix_root)[:, :max_loci], dtype=np.float32)
    gene_cols = np.unique(assignment.edge_col_index)
    gene_cols = gene_cols[gene_cols < max_loci]
    return {
        "assignment": assignment,
        "betas": betas,
        "train_idx": train_idx,
        "test_idx": test_idx,
        "arrays": arrays,
        "sample_ids": sample_ids,
        "class_names": list(class_names),
        "gene_cols": gene_cols,
        "gene_ids": list(assignment.gene_ids),
        "rho_hidden_dim": int(
            (cfg.get("model") or {}).get("rho_hidden_dimension", DEFAULT_RHO_HIDDEN)
        ),
    }


def _gene_indices_for_trait(artifacts: Any, trait: str, gene_ids: list[str]) -> list[int]:
    genes = artifacts.genes
    trait_genes = genes.loc[genes["trait"].astype(str) == trait, "gene_id"].astype(str).tolist()
    id_to_idx = {gid: i for i, gid in enumerate(gene_ids)}
    return [id_to_idx[g] for g in trait_genes if g in id_to_idx]


def _train_kwargs(pack: dict[str, Any], *, arm: str, epochs: int, seed: int, device: str, out: Path, masks: dict[str, Any] | None = None) -> dict[str, Any]:
    kw: dict[str, Any] = dict(
        assignment=pack["assignment"],
        betas=pack["betas"],
        train_idx=pack["train_idx"],
        test_idx=pack["test_idx"],
        ages=pack["arrays"]["age"],
        tissue=pack["arrays"]["tissue"],
        sex=pack["arrays"]["sex"],
        study_ids=pack["arrays"]["study_ids"],
        sample_ids=pack["sample_ids"],
        class_names=pack["class_names"],
        out_dir=out,
        max_epochs=epochs,
        seed=seed,
        device_str=device,
        age_mask=pack["arrays"]["age_mask"],
        tissue_mask=pack["arrays"]["tissue_mask"],
        sex_mask=pack["arrays"]["sex_mask"],
        panel_cols=pack["gene_cols"],
        pool="mean",
        arm=arm,
        rho_hidden_dim=int(pack.get("rho_hidden_dim", DEFAULT_RHO_HIDDEN)),
    )
    if masks:
        kw.update(masks)
    return kw


def run_seed_mask_smoke(*, epochs: int, device: str, max_loci: int) -> dict[str, Any]:
    if not (SEED_PANEL_DIR / "seed_panel.json").is_file():
        raise SystemExit(
            f"missing 9c seed panel at {SEED_PANEL_DIR}; generate via "
            "scripts/run_7g_prime_seed_mask.py --reuse-panels first"
        )
    pack = _load_ats_fold0(max_loci=max_loci)
    artifacts = load_seed_panel(SEED_PANEL_DIR)
    n_genes = len(pack["gene_ids"])
    n_classes = max(len(pack["class_names"]), 2)
    age_idx = _gene_indices_for_trait(artifacts, "age", pack["gene_ids"])
    tissue_idx = _gene_indices_for_trait(artifacts, "tissue", pack["gene_ids"])
    sex_idx = _gene_indices_for_trait(artifacts, "sex", pack["gene_ids"])
    masks = {
        "age_seed_mask": gene_mask_tensor(age_idx, n_genes, n_outputs=1),
        "tissue_seed_mask": gene_mask_tensor(tissue_idx, n_genes, n_outputs=n_classes),
        "sex_seed_mask": gene_mask_tensor(sex_idx, n_genes, n_outputs=2),
    }
    results: dict[str, Any] = {
        "n_seed_genes_age": len(age_idx),
        "n_seed_genes_tissue": len(tissue_idx),
        "n_seed_genes_sex": len(sex_idx),
        "rho_hidden_dim": pack.get("rho_hidden_dim", DEFAULT_RHO_HIDDEN),
    }
    for arm, use_masks in (("G0_onehop", False), ("G1_onehop", True)):
        out = REPORT / "seed_mask" / arm
        payload = train_flat_region_on_arrays(
            **_train_kwargs(
                pack,
                arm=arm,
                epochs=epochs,
                seed=42,
                device=device,
                out=out,
                masks=(masks if use_masks else None),
            )
        )
        m = (payload.get("evaluations") or {}).get("mbs_e2e", {}).get("metrics") or {}
        # flat_region_loop may also put metrics at top level
        if not m and "tissue" in (payload.get("metrics") or {}):
            m = payload["metrics"]
        if not m:
            m = {
                "tissue": {"macro_f1": payload.get("tissue_f1") or payload.get("macro_f1")},
                "age": {"mae": payload.get("age_mae") or payload.get("mae")},
                "sex": {"auroc": payload.get("sex_auroc") or payload.get("auroc")},
            }
        results[arm] = {
            "tissue_f1": (m.get("tissue") or {}).get("macro_f1"),
            "age_mae": (m.get("age") or {}).get("mae"),
            "sex_auroc": (m.get("sex") or {}).get("auroc"),
        }
        print(f"[onehop-smoke] {arm} {results[arm]}", flush=True)
    return results


def run_multiseed_smoke(*, epochs: int, device: str, max_loci: int, seeds: list[int]) -> dict[str, Any]:
    pack = _load_ats_fold0(max_loci=max_loci)
    results: dict[str, Any] = {
        "rho_hidden_dim": pack.get("rho_hidden_dim", DEFAULT_RHO_HIDDEN),
    }
    for seed in seeds:
        arm = f"onehop_mean_s{seed}"
        payload = train_flat_region_on_arrays(
            **_train_kwargs(
                pack,
                arm=arm,
                epochs=epochs,
                seed=seed,
                device=device,
                out=REPORT / "multiseed" / arm,
            )
        )
        m = (payload.get("evaluations") or {}).get("mbs_e2e", {}).get("metrics") or {}
        if not m and "tissue" in (payload.get("metrics") or {}):
            m = payload["metrics"]
        results[arm] = {
            "seed": seed,
            "tissue_f1": (m.get("tissue") or {}).get("macro_f1"),
            "age_mae": (m.get("age") or {}).get("mae"),
            "sex_auroc": (m.get("sex") or {}).get("auroc"),
        }
        print(f"[onehop-smoke] {arm} {results[arm]}", flush=True)
    tissues = [float(v["tissue_f1"]) for v in results.values() if isinstance(v, dict) and v.get("tissue_f1") is not None]
    results["_diversity"] = {
        "tissue_f1_span": (max(tissues) - min(tissues)) if len(tissues) >= 2 else None,
        "n_seeds": len(seeds),
        "ok_distinct": bool(len(tissues) >= 2 and abs(max(tissues) - min(tissues)) > 1e-6),
    }
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--epochs", type=int, default=5)
    parser.add_argument("--max-loci", type=int, default=65536)
    parser.add_argument("--skip-seed-mask", action="store_true")
    parser.add_argument("--skip-multiseed", action="store_true")
    parser.add_argument("--seeds", default="42,43,44")
    args = parser.parse_args()
    REPORT.mkdir(parents=True, exist_ok=True)
    summary: dict[str, Any] = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "epochs": args.epochs,
        "rho_hidden_dim_default": DEFAULT_RHO_HIDDEN,
        "note": (
            "Cheap correctness only. One-hop+seed-mask and one-hop multi-seed "
            "were previously untested; do not scale before these land. "
            f"Uses N-light default rho_hidden={DEFAULT_RHO_HIDDEN}."
        ),
    }
    if not args.skip_seed_mask:
        summary["seed_mask"] = run_seed_mask_smoke(
            epochs=args.epochs, device=args.device, max_loci=args.max_loci
        )
    if not args.skip_multiseed:
        seeds = [int(x) for x in str(args.seeds).split(",") if x.strip()]
        summary["multiseed"] = run_multiseed_smoke(
            epochs=args.epochs, device=args.device, max_loci=args.max_loci, seeds=seeds
        )
    (REPORT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    lines = [
        "# One-hop correctness smokes (Milestone 10)",
        "",
        f"Generated: `{summary['generated_at']}`",
        "",
        summary["note"],
        "",
    ]
    if "seed_mask" in summary:
        lines += ["## Seed-mask (fold 0, G0 vs G1 on one-hop)", ""]
        for arm, row in summary["seed_mask"].items():
            lines.append(f"- **{arm}**: `{row}`")
        lines.append("")
        g0 = summary["seed_mask"].get("G0_onehop") or {}
        g1 = summary["seed_mask"].get("G1_onehop") or {}
        lines += [
            "### Interpretation",
            "",
            "G0 = dense one-hop (no seed mask); G1 = trait-associated seed genes only. "
            "On cascade (9c), G0 beat G1–G3 for age-primary training. Here we ask whether "
            "one-hop’s single-stage pooling interacts differently with sparse masking. "
            "Stable train + sensible metrics = correctness pass; adopt masking only if G1 "
            "clearly helps on the primary tasks.",
            "",
            f"- G0 tissue/age/sex: `{g0}`",
            f"- G1 tissue/age/sex: `{g1}`",
            "",
        ]
    if "multiseed" in summary:
        lines += ["## Multi-seed restarts (fold 0)", ""]
        for arm, row in summary["multiseed"].items():
            lines.append(f"- **{arm}**: `{row}`")
        lines.append("")
        div = summary["multiseed"].get("_diversity") or {}
        lines += [
            "### Interpretation",
            "",
            "After the seed-offset fix, restarts must produce **distinct** sensible runs "
            f"(tissue F1 span `{div.get('tissue_f1_span')}`, ok_distinct=`{div.get('ok_distinct')}`). "
            "A zero span means the seed path is still broken; do not trust Milestone 12 "
            "restart ensembling until this is true.",
            "",
        ]
    (REPORT / "analysis.md").write_text("\n".join(lines) + "\n")
    print(f"wrote {REPORT / 'analysis.md'}", flush=True)


if __name__ == "__main__":
    main()
