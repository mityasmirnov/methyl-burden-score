#!/usr/bin/env python3
"""Build and freeze a study-grouped k-fold split for the nine-pack Hub cohort.

Mirrors the mechanism already used for ``hub-ats-7e-3fold-v1``
(``dev_cv.samples_from_phenotype_table`` + ``dev_cv.freeze_outer_folds``,
``evaluation.splits.build_outer_study_grouped_folds``) but points at the
nine-pack phenotype table / tissue ontology (34,234 samples: age, tissue,
sex, plus disease/cancer/blood/brain mask flags) instead of the ATS-only
table. No new split logic — this is the existing, leakage-checked
(``assert_no_study_leakage``) builder, applied to a bigger cohort.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from mbs.paths import DataPaths
from mbs.training.dev_cv import freeze_outer_folds, samples_from_phenotype_table

DEFAULT_SPLIT_ID = "hub-nine-pack-3fold-v1"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--phenotype-table",
        default="canonical/phenotypes/sample_phenotype_table_hub_nine_pack_v1.parquet",
    )
    parser.add_argument(
        "--tissue-ontology",
        default="canonical/phenotypes/tissue_ontology_hub_nine_pack_v1.yaml",
    )
    parser.add_argument("--split-id", default=DEFAULT_SPLIT_ID)
    parser.add_argument("--n-folds", type=int, default=3)
    parser.add_argument("--val-fraction", type=float, default=0.15)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    paths = DataPaths.from_environment()
    table_path = paths.data_root / args.phenotype_table
    ontology_path = paths.data_root / args.tissue_ontology
    out_dir = paths.artifact_root / "splits" / args.split_id

    print(f"[split] loading samples from {table_path}", flush=True)
    samples, phenotypes = samples_from_phenotype_table(table_path, ontology_path=ontology_path)
    print(f"[split] {len(samples)} samples, {len({s['study_id'] for s in samples})} studies", flush=True)

    if out_dir.exists():
        raise SystemExit(
            f"{out_dir} already exists -- refusing to overwrite a frozen split. "
            "Delete it first if you really mean to rebuild."
        )

    fold_pack = freeze_outer_folds(
        samples,
        out_dir=out_dir,
        n_folds=args.n_folds,
        seed=args.seed,
        val_fraction=args.val_fraction,
        split_id=args.split_id,
    )

    print(f"[split] wrote {out_dir / 'folds.json'} (sha256={fold_pack['sha256'][:16]}…)", flush=True)
    for i, fold in enumerate(fold_pack["folds"]):
        n_train = len(fold["train_sample_ids"])
        n_val = len(fold["validation_sample_ids"])
        n_test = len(fold.get("external_test_sample_ids") or [])
        print(f"[split] fold {i}: train={n_train} val={n_val} test={n_test}", flush=True)

    # Per-trait mask coverage, so we know which traits are actually usable
    # per fold before wiring any training config to this split.
    by_id = {p.sample_id: p for p in phenotypes}
    for i, fold in enumerate(fold_pack["folds"]):
        test_ids = fold.get("external_test_sample_ids") or []
        masks = ["age_mask", "tissue_mask", "sex_mask"]
        counts = {m: sum(1 for sid in test_ids if getattr(by_id.get(sid), m, False)) for m in masks}
        print(f"[split] fold {i} external_test mask coverage: {counts}", flush=True)

    print("[split] done", flush=True)


if __name__ == "__main__":
    main()
