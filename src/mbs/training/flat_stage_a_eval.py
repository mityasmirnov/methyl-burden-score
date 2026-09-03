"""Stage A evaluation suite for FlatDeepSetRegion (N-light-gene-*) arms.

Inline screen (P2/P4-like): ``mbs_e2e`` + ``mbs_linear_probe`` (test split only).
``mbs_enet`` is **post-hoc** via ``scripts/eval_mbs_enet_from_scores.py`` so the
GPU queue is not blocked on sklearn elastic-net. One-hop models have no region
RBS layer, so ``rbs_*`` / orphan fusion modes are omitted.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import torch
from torch import nn

from mbs.models import FlatDeepSet
from mbs.training.cascade_scores import fusion_feature_matrix
from mbs.training.dataset import FlatSampleRecord, pack_records_to_batch
from mbs.training.late_fusion import evaluate_late_fusion
from mbs.training.multitask import MultitaskHeads
from mbs.training.phenotypes import SamplePhenotype
from mbs.training.transparent_baselines import (
    evaluate_multitask_predictions,
    run_elasticnet_multitask,
)


def score_flat_mbs_matrix(
    *,
    phenotypes: list[SamplePhenotype],
    materialize_fn,
    model: FlatDeepSet,
    device: torch.device,
    n_genes: int,
    batch_size: int = 16,
) -> tuple[np.ndarray, np.ndarray]:
    """Encode phenotypes → MBS/present matrices with GPU micro-batches."""
    if not phenotypes:
        return (
            np.zeros((0, max(n_genes, 1)), dtype=np.float32),
            np.zeros((0, max(n_genes, 1)), dtype=bool),
        )
    model.eval()
    mbs_chunks: list[np.ndarray] = []
    present_chunks: list[np.ndarray] = []
    step = max(1, int(batch_size))
    for start in range(0, len(phenotypes), step):
        chunk_ph = phenotypes[start : start + step]
        chunk: list[FlatSampleRecord] = [materialize_fn(ph) for ph in chunk_ph]
        batch = pack_records_to_batch(
            chunk,
            n_genes=n_genes,
            age_values=[None] * len(chunk),
            age_enabled=[False] * len(chunk),
            tissue_enabled=[False] * len(chunk),
            sex_enabled=[False] * len(chunk),
            sex_class_indices=[0] * len(chunk),
        ).to(device, non_blocking=device.type == "cuda")
        with torch.no_grad():
            n_instances = len(batch.sample_ids) * int(batch.n_genes)
            out = model(batch.cpg_features, batch.cpg_to_gene, n_instances)
            mbs = out["mbs"].view(len(batch.sample_ids), batch.n_genes)
            present = out["present"].view(len(batch.sample_ids), batch.n_genes)
        mbs_chunks.append(mbs.detach().float().cpu().numpy())
        present_chunks.append(present.detach().cpu().numpy())
    return np.concatenate(mbs_chunks, axis=0), np.concatenate(present_chunks, axis=0)


def _phenotype_arrays(
    phenotypes: list[SamplePhenotype],
) -> dict[str, np.ndarray]:
    return {
        "age": np.asarray([float(p.age or 0.0) for p in phenotypes], dtype=np.float64),
        "age_mask": np.asarray([bool(p.age_mask) for p in phenotypes], dtype=bool),
        "tissue": np.asarray([int(p.class_index) for p in phenotypes], dtype=np.int64),
        "tissue_mask": np.asarray([bool(p.tissue_mask) for p in phenotypes], dtype=bool),
        "sex": np.asarray([int(p.sex_class_index or 0) for p in phenotypes], dtype=np.int64),
        "sex_mask": np.asarray([bool(p.sex_mask) for p in phenotypes], dtype=bool),
        "study_ids": np.asarray([str(p.study_id or p.sample_id) for p in phenotypes], dtype=object),
    }


def evaluate_flat_mbs_e2e(
    *,
    heads: MultitaskHeads,
    mbs_test: np.ndarray,
    present_test: np.ndarray,
    phenotypes_test: list[SamplePhenotype],
    phenotypes_train: list[SamplePhenotype],
    class_names: list[str],
    device: torch.device,
    age_mean: float,
    age_std: float,
) -> dict[str, Any]:
    """End-to-end MultitaskHeads on test MBS only (destandardize age to years)."""
    if not phenotypes_test:
        raise ValueError("evaluate_flat_mbs_e2e requires test phenotypes")
    heads.eval()
    mbs_t = torch.from_numpy(np.asarray(mbs_test, dtype=np.float32)).to(device)
    present_t = torch.from_numpy(np.asarray(present_test, dtype=bool)).to(device)
    with torch.no_grad():
        age_hat = heads.forward_age(mbs_t, present_t).detach().cpu().numpy().reshape(-1)
        age_hat = age_hat * float(age_std) + float(age_mean)
        tissue_logits = heads.forward_tissue(mbs_t, present_t)
        tissue_pred = tissue_logits.detach().cpu().numpy().argmax(axis=1)
        tissue_proba = torch.softmax(tissue_logits, dim=-1).detach().cpu().numpy()
        sex_pred: np.ndarray | None = None
        sex_proba: np.ndarray | None = None
        if heads.sex_head is not None:
            sex_logits = heads.forward_sex(mbs_t, present_t)
            if sex_logits is not None:
                sex_pred = sex_logits.detach().cpu().numpy().argmax(axis=1)
                sex_proba = torch.softmax(sex_logits, dim=-1).detach().cpu().numpy()
    preds: dict[str, np.ndarray | None] = {
        "age": age_hat,
        "tissue": tissue_pred,
        "tissue_proba": tissue_proba,
        "tissue_classes": np.arange(tissue_proba.shape[1], dtype=np.int64),
        "sex": sex_pred,
    }
    if sex_proba is not None:
        preds["sex_proba"] = sex_proba
        preds["sex_classes"] = np.arange(sex_proba.shape[1], dtype=np.int64)
    train_arr = _phenotype_arrays(phenotypes_train)
    test_arr = _phenotype_arrays(phenotypes_test)
    tissue_valid_classes = None
    if train_arr["tissue_mask"].any():
        tissue_valid_classes = set(
            train_arr["tissue"][train_arr["tissue_mask"]].astype(int).tolist()
        )
    metrics = evaluate_multitask_predictions(
        preds=preds,
        age=test_arr["age"],
        age_mask=test_arr["age_mask"],
        tissue=test_arr["tissue"],
        tissue_mask=test_arr["tissue_mask"],
        sex=test_arr["sex"],
        sex_mask=test_arr["sex_mask"],
        study_ids=test_arr["study_ids"],
        tissue_class_names=list(class_names) if class_names else None,
        tissue_valid_classes=tissue_valid_classes,
    )
    return {
        "metrics": metrics,
        "evaluation": "mbs_e2e",
        "eval_split": "test",
        "n_eval_samples": int(len(phenotypes_test)),
        "n_score_features": int(mbs_test.shape[1]),
    }


def evaluate_flat_mbs_linear_probe(
    *,
    mbs_all: np.ndarray,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    arrays: dict[str, np.ndarray],
    class_names: list[str],
) -> dict[str, Any]:
    """CPU linear probe on frozen MBS (same contract as cascade ``mbs_linear_probe``)."""
    blocks = {"mbs": np.asarray(mbs_all, dtype=np.float32)}
    x = fusion_feature_matrix(blocks, mode="mbs_only")
    out = evaluate_late_fusion(
        scores_train=x[train_idx],
        scores_test=x[test_idx],
        age_train=arrays["age"][train_idx],
        age_mask_train=arrays["age_mask"][train_idx],
        tissue_train=arrays["tissue"][train_idx],
        tissue_mask_train=arrays["tissue_mask"][train_idx],
        sex_train=arrays["sex"][train_idx],
        sex_mask_train=arrays["sex_mask"][train_idx],
        age_test=arrays["age"][test_idx],
        age_mask_test=arrays["age_mask"][test_idx],
        tissue_test=arrays["tissue"][test_idx],
        tissue_mask_test=arrays["tissue_mask"][test_idx],
        sex_test=arrays["sex"][test_idx],
        sex_mask_test=arrays["sex_mask"][test_idx],
        study_ids_test=arrays["study_ids"][test_idx],
        tissue_class_names=list(class_names) if class_names else None,
        # Match cascade P2/P4 linear probe (default lbfgs LogReg).
        fusion=None,
    )
    out["evaluation"] = "mbs_linear_probe"
    out["eval_split"] = "test"
    out["fusion_block_mode"] = "mbs_only"
    out["n_eval_samples"] = int(np.asarray(test_idx).size)
    return out


def evaluate_flat_mbs_enet(
    *,
    mbs_all: np.ndarray,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    arrays: dict[str, np.ndarray],
    class_names: list[str],
) -> dict[str, Any]:
    """CPU elastic-net readout on frozen MBS (same contract as cascade ``mbs_enet``)."""
    blocks = {"mbs": np.asarray(mbs_all, dtype=np.float32)}
    x = fusion_feature_matrix(blocks, mode="mbs_only")
    out = run_elasticnet_multitask(
        x_train=x[train_idx],
        x_test=x[test_idx],
        age_train=arrays["age"][train_idx],
        age_mask_train=arrays["age_mask"][train_idx],
        tissue_train=arrays["tissue"][train_idx],
        tissue_mask_train=arrays["tissue_mask"][train_idx],
        sex_train=arrays["sex"][train_idx],
        sex_mask_train=arrays["sex_mask"][train_idx],
        age_test=arrays["age"][test_idx],
        age_mask_test=arrays["age_mask"][test_idx],
        tissue_test=arrays["tissue"][test_idx],
        tissue_mask_test=arrays["tissue_mask"][test_idx],
        sex_test=arrays["sex"][test_idx],
        sex_mask_test=arrays["sex_mask"][test_idx],
        study_ids_test=arrays["study_ids"][test_idx],
        tissue_class_names=list(class_names) if class_names else None,
    )
    out["evaluation"] = "mbs_enet"
    out["eval_split"] = "test"
    out["n_eval_samples"] = int(np.asarray(test_idx).size)
    out["n_score_features"] = int(out.get("n_features", x.shape[1]))
    return out


def _run_inline_cpu_probes(
    *,
    mbs_all: np.ndarray,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    arrays: dict[str, np.ndarray],
    class_names: list[str],
    include_mbs_enet: bool = False,
) -> dict[str, Any]:
    """Fit Stage A CPU probes. Default: linear only (enet is post-hoc)."""
    print("[flat] Stage A CPU probe: mbs_linear_probe…", flush=True)  # noqa: T201
    out: dict[str, Any] = {
        "mbs_linear_probe": evaluate_flat_mbs_linear_probe(
            mbs_all=mbs_all,
            train_idx=train_idx,
            test_idx=test_idx,
            arrays=arrays,
            class_names=list(class_names),
        )
    }
    if include_mbs_enet:
        print("[flat] Stage A CPU probe: mbs_enet (inline; prefer post-hoc)…", flush=True)  # noqa: T201
        out["mbs_enet"] = evaluate_flat_mbs_enet(
            mbs_all=mbs_all,
            train_idx=train_idx,
            test_idx=test_idx,
            arrays=arrays,
            class_names=list(class_names),
        )
    return out


def save_stage_a_probe_inputs(
    score_dir: Path,
    *,
    mbs_all: np.ndarray,
    present_all: np.ndarray,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    arrays: dict[str, np.ndarray],
    class_names: list[str],
    sample_ids: list[str],
) -> None:
    """Persist frozen-MBS probe inputs so CPU eval can overlap the next GPU fold."""
    score_dir.mkdir(parents=True, exist_ok=True)
    np.save(score_dir / "mbs.npy", mbs_all)
    np.save(score_dir / "mbs_present.npy", present_all.astype(np.uint8))
    np.savez_compressed(
        score_dir / "stage_a_probe_inputs.npz",
        train_idx=np.asarray(train_idx, dtype=np.int64),
        test_idx=np.asarray(test_idx, dtype=np.int64),
        age=np.asarray(arrays["age"], dtype=np.float64),
        age_mask=np.asarray(arrays["age_mask"], dtype=bool),
        tissue=np.asarray(arrays["tissue"], dtype=np.int64),
        tissue_mask=np.asarray(arrays["tissue_mask"], dtype=bool),
        sex=np.asarray(arrays["sex"], dtype=np.int64),
        sex_mask=np.asarray(arrays["sex_mask"], dtype=bool),
        study_ids=np.asarray(arrays["study_ids"], dtype=object),
    )
    import json

    (score_dir / "sample_ids.json").write_text(
        json.dumps(list(sample_ids), indent=2) + "\n", encoding="utf-8"
    )
    (score_dir / "stage_a_probe_meta.json").write_text(
        json.dumps({"class_names": list(class_names), "pending_cpu_probes": True}, indent=2)
        + "\n",
        encoding="utf-8",
    )


def complete_flat_stage_a_cpu_probes(
    run_root: str | Path,
    *,
    include_mbs_enet: bool = False,
) -> dict[str, Any]:
    """Finish deferred ``mbs_linear_probe`` (and optional enet) into metrics.json."""
    import json

    run_root_p = Path(run_root)
    score_dir = run_root_p / "scores"
    metrics_path = run_root_p / "metrics.json"
    inputs_path = score_dir / "stage_a_probe_inputs.npz"
    meta_path = score_dir / "stage_a_probe_meta.json"
    mbs_path = score_dir / "mbs.npy"
    if not inputs_path.is_file() or not mbs_path.is_file():
        raise FileNotFoundError(f"missing Stage A probe inputs under {score_dir}")
    blob = np.load(inputs_path, allow_pickle=True)
    mbs_all = np.load(mbs_path)
    arrays = {
        "age": blob["age"],
        "age_mask": blob["age_mask"],
        "tissue": blob["tissue"],
        "tissue_mask": blob["tissue_mask"],
        "sex": blob["sex"],
        "sex_mask": blob["sex_mask"],
        "study_ids": blob["study_ids"],
    }
    class_names = ["tissue"]
    if meta_path.is_file():
        class_names = list(json.loads(meta_path.read_text(encoding="utf-8")).get("class_names") or class_names)
    metrics = {}
    if metrics_path.is_file():
        metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
    if "include_mbs_enet" in metrics:
        include_mbs_enet = bool(metrics.get("include_mbs_enet"))
    probes = _run_inline_cpu_probes(
        mbs_all=mbs_all,
        train_idx=blob["train_idx"],
        test_idx=blob["test_idx"],
        arrays=arrays,
        class_names=class_names,
        include_mbs_enet=include_mbs_enet,
    )
    evaluations = dict(metrics.get("evaluations") or {})
    evaluations.update(probes)
    metrics["evaluations"] = evaluations
    metrics["cpu_probes_pending"] = False
    metrics["include_mbs_enet"] = bool(include_mbs_enet)
    metrics_path.write_text(json.dumps(metrics, indent=2, default=str) + "\n", encoding="utf-8")
    if meta_path.is_file():
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        meta["pending_cpu_probes"] = False
        meta_path.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    print(  # noqa: T201
        f"[flat] deferred CPU probes complete for {run_root_p.name}: "
        + ", ".join(sorted(probes.keys())),
        flush=True,
    )
    return probes


def build_stage_a_flat_evaluations(
    *,
    model: FlatDeepSet,
    head: nn.Module,
    train_phenotypes: list[SamplePhenotype],
    val_phenotypes: list[SamplePhenotype] | None,
    test_phenotypes: list[SamplePhenotype],
    materialize_fn,
    device: torch.device,
    n_genes: int,
    class_names: list[str],
    age_mean: float,
    age_std: float,
    batch_size: int,
    score_dir: Path | None = None,
    defer_cpu_probes: bool = False,
    include_mbs_enet: bool = False,
) -> dict[str, Any]:
    """Score train/val/test MBS and return Stage A evaluation dict.

    Default inline suite matches P2/P4 screen: ``mbs_e2e`` + ``mbs_linear_probe``.
    ``mbs_enet`` stays off unless ``include_mbs_enet`` (prefer post-hoc script).
    When ``defer_cpu_probes`` is true, only ``mbs_e2e`` runs here; linear inputs
    are written under ``score_dir`` for ``complete_flat_stage_a_cpu_probes``.
    """
    if not isinstance(head, MultitaskHeads):
        raise TypeError("Stage A flat evaluations require MultitaskHeads")
    val_ph = list(val_phenotypes or [])
    ordered = list(train_phenotypes) + val_ph + list(test_phenotypes)
    print(  # noqa: T201
        f"[flat] Stage A scoring MBS n={len(ordered)} "
        f"(train={len(train_phenotypes)} val={len(val_ph)} test={len(test_phenotypes)}) "
        f"batch_size={batch_size}",
        flush=True,
    )
    mbs_all, present_all = score_flat_mbs_matrix(
        phenotypes=ordered,
        materialize_fn=materialize_fn,
        model=model,
        device=device,
        n_genes=n_genes,
        batch_size=batch_size,
    )
    n_train = len(train_phenotypes)
    n_val = len(val_ph)
    train_idx = np.arange(0, n_train, dtype=np.int64)
    test_idx = np.arange(n_train + n_val, n_train + n_val + len(test_phenotypes), dtype=np.int64)
    arrays = _phenotype_arrays(ordered)
    mbs_te = mbs_all[test_idx]
    present_te = present_all[test_idx]
    if score_dir is not None:
        save_stage_a_probe_inputs(
            score_dir,
            mbs_all=mbs_all,
            present_all=present_all,
            train_idx=train_idx,
            test_idx=test_idx,
            arrays=arrays,
            class_names=list(class_names),
            sample_ids=[p.sample_id for p in ordered],
        )
    evaluations: dict[str, Any] = {
        "mbs_e2e": evaluate_flat_mbs_e2e(
            heads=head,
            mbs_test=mbs_te,
            present_test=present_te,
            phenotypes_test=list(test_phenotypes),
            phenotypes_train=list(train_phenotypes),
            class_names=list(class_names),
            device=device,
            age_mean=age_mean,
            age_std=age_std,
        ),
    }
    if defer_cpu_probes:
        print(  # noqa: T201
            "[flat] deferring CPU linear probe to overlap next fold GPU train "
            f"(mbs_enet inline={'on' if include_mbs_enet else 'off → post-hoc'})",
            flush=True,
        )
        return {
            "evaluations": evaluations,
            "primary_evaluation": "mbs_e2e",
            "gene_linked_only": True,
            "n_genes": int(n_genes),
            "score_dir": None if score_dir is None else str(score_dir),
            "cpu_probes_pending": True,
            "include_mbs_enet": bool(include_mbs_enet),
            "eval_split": "test",
        }
    evaluations.update(
        _run_inline_cpu_probes(
            mbs_all=mbs_all,
            train_idx=train_idx,
            test_idx=test_idx,
            arrays=arrays,
            class_names=list(class_names),
            include_mbs_enet=include_mbs_enet,
        )
    )
    print(  # noqa: T201
        "[flat] Stage A evaluations ready: "
        + ", ".join(sorted(evaluations.keys())),
        flush=True,
    )
    return {
        "evaluations": evaluations,
        "primary_evaluation": "mbs_e2e",
        "gene_linked_only": True,
        "n_genes": int(n_genes),
        "score_dir": None if score_dir is None else str(score_dir),
        "cpu_probes_pending": False,
        "include_mbs_enet": bool(include_mbs_enet),
        "eval_split": "test",
    }
