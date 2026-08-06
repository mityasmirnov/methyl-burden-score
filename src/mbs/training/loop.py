"""Plain PyTorch training loop for the flat DeepRVAT-style baseline."""

from __future__ import annotations

import json
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import torch
import torch.nn.functional as F
import yaml
from torch import nn

from mbs.evaluation.splits import build_study_grouped_split, partition_studies_by_sample_count
from mbs.matrix.store import (
    matrix_store_paths,
    open_betas_zarr,
    read_locus_index,
    read_sample_index,
)
from mbs.models import FlatDeepSet, SeedMaskedLinearHead
from mbs.segment_ops import PoolName
from mbs.static_features.store import (
    open_embeddings_zarr,
    read_loci_index,
    static_feature_store_paths,
)
from mbs.training.dataset import (
    FlatBatch,
    FlatSampleRecord,
    build_flat_sample,
    make_synthetic_overfit_bundle,
    make_synthetic_study_holdout_bundle,
    pack_records_to_batch,
    record_to_batch,
)
from mbs.training.features import build_static_column_table
from mbs.training.locus_gene import LocusGeneIndex, build_locus_gene_index, load_graph_tables
from mbs.training.multitask import MultitaskHeads, masked_multitask_loss
from mbs.training.phenotype_table import load_tissue_ontology
from mbs.training.phenotypes import (
    SamplePhenotype,
    load_gse35069_phenotypes,
    load_hub_regression_phenotypes,
    load_hub_sample_info_phenotypes,
    load_multitask_phenotypes,
)
from mbs.training.run_artifacts import (
    checkpoint_dir,
    collect_environment,
    config_sha256,
    run_dir,
    save_checkpoint,
    write_run_artifacts,
)


@dataclass(frozen=True, slots=True)
class TrainResult:
    run_id: str
    run_dir: Path
    checkpoint_dir: Path
    metrics: dict[str, Any]
    best_epoch: int
    tensorboard_url: str | None = None
    tensorboard_port: int | None = None
    monitor_hint: str | None = None


@dataclass(slots=True)
class _PilotStore:
    phenotypes: list[SamplePhenotype]
    sample_row_by_id: dict[str, int]
    betas: Any
    static_by_col: np.ndarray
    static_valid: np.ndarray
    locus_gene: LocusGeneIndex
    epsilon: float
    n_cols: int


def resolve_device(device_str: str, *, require_cuda: bool = False) -> torch.device:
    requested = device_str.strip().lower()
    if requested.startswith("cuda"):
        if not torch.cuda.is_available():
            raise RuntimeError(
                "CUDA was requested but torch.cuda.is_available() is False. "
                "Check the driver / cu128 torch wheel and CUDA_VISIBLE_DEVICES."
            )
        if torch.cuda.device_count() < 1:
            raise RuntimeError("CUDA requested but no visible devices (check CUDA_VISIBLE_DEVICES)")
        return torch.device("cuda:0")
    if require_cuda:
        raise RuntimeError("require_cuda=True but device is not cuda")
    if requested == "cpu":
        return torch.device("cpu")
    return torch.device(device_str)


def load_experiment_config(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        raise TypeError(f"config must be a mapping: {path}")
    return data


def _set_seed(seed: int) -> None:
    torch.manual_seed(seed)
    np.random.seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def _class_weights(labels: list[int], n_classes: int) -> torch.Tensor:
    counts = np.bincount(np.asarray(labels, dtype=np.int64), minlength=n_classes).astype(np.float64)
    counts = np.maximum(counts, 1.0)
    weights = counts.sum() / (n_classes * counts)
    return torch.tensor(weights, dtype=torch.float32)


def _packed_mbs(model: FlatDeepSet, batch: FlatBatch) -> tuple[torch.Tensor, torch.Tensor]:
    """Run FlatDeepSet on a (possibly multi-sample) packed batch → ``[B, G]``."""
    batch_size = len(batch.sample_ids)
    n_instances = batch_size * int(batch.n_genes)
    output = model(batch.cpg_features, batch.cpg_to_gene, n_instances)
    mbs = output["mbs"].view(batch_size, batch.n_genes)
    present = output["present"].view(batch_size, batch.n_genes)
    return mbs, present


def _forward_tissue_loss(
    *,
    model: FlatDeepSet,
    head: SeedMaskedLinearHead,
    batch: FlatBatch,
    class_weights: torch.Tensor | None,
) -> tuple[torch.Tensor, dict[str, float]]:
    mbs, present = _packed_mbs(model, batch)
    logits = head(mbs, present)
    targets = batch.tissue_target.reshape(-1)
    tissue_loss = F.cross_entropy(
        logits,
        targets,
        weight=class_weights.to(logits.device) if class_weights is not None else None,
    )
    pred = logits.argmax(dim=-1)
    metrics = {
        "loss": float(tissue_loss.detach().item()),
        "tissue_correct": float((pred == targets).sum().item()),
        "tissue_n": float(targets.numel()),
        "mae": 0.0,
        "age_n": 0.0,
    }
    return tissue_loss, metrics


def _forward_age_loss(
    *,
    model: FlatDeepSet,
    head: nn.Module,
    batch: FlatBatch,
) -> tuple[torch.Tensor, dict[str, float]]:
    if batch.age_target is None:
        raise RuntimeError("age_target required for regression task")
    mbs, present = _packed_mbs(model, batch)
    pred = head(mbs * present.to(dtype=mbs.dtype)).squeeze(-1)
    target = batch.age_target.reshape(pred.shape)
    loss = F.mse_loss(pred, target)
    mae = float((pred.detach() - target.detach()).abs().mean().item())
    return loss, {
        "loss": float(loss.detach().item()),
        "tissue_correct": 0.0,
        "tissue_n": 0.0,
        "mae": mae,
        "age_n": float(pred.numel()),
    }


def _materialize_record(
    phenotype: SamplePhenotype,
    store: _PilotStore,
) -> FlatSampleRecord:
    row = store.sample_row_by_id[phenotype.sample_id]
    beta_row = np.asarray(store.betas[row, : store.n_cols], dtype=np.float32)
    return build_flat_sample(
        phenotype=phenotype,
        beta_row=beta_row,
        static_by_col=store.static_by_col,
        static_valid=store.static_valid,
        locus_gene=store.locus_gene,
        epsilon=store.epsilon,
    )


def _label_flags_for_record(
    *,
    record: FlatSampleRecord,
    ph: SamplePhenotype | None,
    task: str,
    age_mean: float,
    age_std: float,
) -> tuple[float | None, bool, bool, bool, int]:
    """Return age value, age/tissue/sex enabled flags, and sex class index."""
    if task == "multitask":
        age_enabled = bool(ph.age_mask) if ph is not None else False
        tissue_enabled = bool(ph.tissue_mask) if ph is not None else False
        sex_enabled = bool(ph.sex_mask) if ph is not None else False
        sex_cls = int(ph.sex_class_index) if ph is not None else 0
        age_value = None
        if age_enabled:
            if ph is None or ph.age is None:
                raise RuntimeError(f"missing age for sample {record.sample_id}")
            age_value = (float(ph.age) - age_mean) / age_std
        return age_value, age_enabled, tissue_enabled, sex_enabled, sex_cls
    if task == "regression":
        if ph is None or ph.age is None:
            raise RuntimeError(f"missing age for sample {record.sample_id}")
        return (float(ph.age) - age_mean) / age_std, True, False, False, 0
    return None, False, True, False, 0


def _run_epoch(
    *,
    records: list[FlatSampleRecord] | None,
    phenotypes: list[SamplePhenotype] | None,
    pilot_store: _PilotStore | None,
    model: FlatDeepSet,
    head: nn.Module,
    optimizer: torch.optim.Optimizer | None,
    device: torch.device,
    n_genes: int,
    class_weights: torch.Tensor | None,
    use_amp: bool,
    grad_clip: float,
    train: bool,
    task: str = "multiclass",
    age_mean: float = 0.0,
    age_std: float = 1.0,
    lambda_age: float = 1.0,
    lambda_tissue: float = 1.0,
    lambda_sex: float = 1.0,
    huber_delta: float = 1.0,
    age_loss_name: str = "huber",
    batch_size: int = 1,
) -> dict[str, float]:
    if train:
        model.train()
        head.train()
    else:
        model.eval()
        head.eval()

    total_loss = 0.0
    total_correct = 0.0
    total_mae_sum = 0.0
    age_n = 0.0
    tissue_n = 0.0
    sex_n = 0.0
    sex_correct = 0.0
    n = 0
    amp_dtype = torch.bfloat16 if device.type == "cuda" else torch.float32
    context = torch.enable_grad() if train else torch.no_grad()
    step = max(1, int(batch_size))

    ph_by_id = {p.sample_id: p for p in (phenotypes or [])}
    if records is not None:
        all_records: list[FlatSampleRecord] | None = records
        all_phenotypes: list[SamplePhenotype] | None = None
    else:
        all_records = None
        if phenotypes is None or pilot_store is None:
            raise RuntimeError("pilot phenotypes and store are required when records is None")
        all_phenotypes = phenotypes
    n_items = len(all_records) if all_records is not None else len(all_phenotypes or [])

    with context:
        for start in range(0, n_items, step):
            if all_records is not None:
                chunk = all_records[start : start + step]
            else:
                if all_phenotypes is None or pilot_store is None:
                    raise RuntimeError("missing phenotypes/store for on-the-fly materialization")
                chunk = [
                    _materialize_record(ph, pilot_store)
                    for ph in all_phenotypes[start : start + step]
                ]
            age_values: list[float | None] = []
            age_flags: list[bool] = []
            tissue_flags: list[bool] = []
            sex_flags: list[bool] = []
            sex_idxs: list[int] = []
            for record in chunk:
                ph = ph_by_id.get(record.sample_id)
                age_value, age_on, tissue_on, sex_on, sex_cls = _label_flags_for_record(
                    record=record,
                    ph=ph,
                    task=task,
                    age_mean=age_mean,
                    age_std=age_std,
                )
                age_values.append(age_value)
                age_flags.append(age_on)
                tissue_flags.append(tissue_on)
                sex_flags.append(sex_on)
                sex_idxs.append(sex_cls)
            if step == 1:
                batch = record_to_batch(
                    chunk[0],
                    n_genes=n_genes,
                    age_value=age_values[0],
                    age_enabled=age_flags[0],
                    tissue_enabled=tissue_flags[0],
                    sex_enabled=sex_flags[0],
                    sex_class_index=sex_idxs[0],
                ).to(device)
            else:
                batch = pack_records_to_batch(
                    chunk,
                    n_genes=n_genes,
                    age_values=age_values,
                    age_enabled=age_flags,
                    tissue_enabled=tissue_flags,
                    sex_enabled=sex_flags,
                    sex_class_indices=sex_idxs,
                ).to(device)
            if train and optimizer is not None:
                optimizer.zero_grad(set_to_none=True)
            with torch.autocast(device_type=device.type, dtype=amp_dtype, enabled=use_amp):
                if task == "multitask":
                    if not isinstance(head, MultitaskHeads):
                        raise TypeError("multitask task requires MultitaskHeads")
                    mbs, present = _packed_mbs(model, batch)
                    result = masked_multitask_loss(
                        mbs=mbs,
                        present=present,
                        heads=head,
                        batch=batch,
                        lambda_age=lambda_age,
                        lambda_tissue=lambda_tissue,
                        lambda_sex=lambda_sex,
                        huber_delta=huber_delta,
                        age_loss=age_loss_name,
                        class_weights=class_weights,
                    )
                    loss = result.loss
                    metrics = result.metrics
                elif task == "regression":
                    loss, metrics = _forward_age_loss(model=model, head=head, batch=batch)
                else:
                    loss, metrics = _forward_tissue_loss(
                        model=model,
                        head=head,  # type: ignore[arg-type]
                        batch=batch,
                        class_weights=class_weights,
                    )
            if train and optimizer is not None:
                loss.backward()
                if grad_clip > 0:
                    nn.utils.clip_grad_norm_(
                        list(model.parameters()) + list(head.parameters()),
                        grad_clip,
                    )
                optimizer.step()
            batch_n = float(len(chunk))
            total_loss += metrics["loss"] * batch_n
            step_age_n = float(metrics.get("age_n", 1.0 if task == "regression" else 0.0))
            step_tissue_n = float(metrics.get("tissue_n", 0.0 if task == "regression" else batch_n))
            total_correct += metrics["tissue_correct"]
            total_mae_sum += metrics["mae"] * step_age_n
            age_n += step_age_n
            tissue_n += step_tissue_n
            sex_n += float(metrics.get("sex_n", 0.0))
            sex_correct += float(metrics.get("sex_correct", 0.0))
            n += int(batch_n)

    return {
        "loss": total_loss / max(n, 1),
        "accuracy": total_correct / max(tissue_n, 1.0),
        "sex_accuracy": sex_correct / max(sex_n, 1.0),
        "mae": total_mae_sum / max(age_n, 1.0),
        "n_samples": float(n),
        "age_n": age_n,
        "tissue_n": tissue_n,
        "sex_n": sex_n,
    }


def _split_by_donor(
    phenotypes: list[SamplePhenotype],
    *,
    train_donors: set[str],
    val_donors: set[str],
) -> tuple[list[SamplePhenotype], list[SamplePhenotype], dict[str, Any]]:
    train_ph: list[SamplePhenotype] = []
    val_ph: list[SamplePhenotype] = []
    for ph in phenotypes:
        if ph.donor_id in train_donors:
            train_ph.append(ph)
        elif ph.donor_id in val_donors:
            val_ph.append(ph)
        else:
            raise ValueError(f"donor {ph.donor_id!r} not in train or val sets")
    split = {
        "train_donors": sorted(train_donors),
        "val_donors": sorted(val_donors),
        "train_sample_ids": [p.sample_id for p in train_ph],
        "val_sample_ids": [p.sample_id for p in val_ph],
    }
    return train_ph, val_ph, split


def train_flat_baseline(
    *,
    project_root: Path,
    data_root: Path,
    artifact_root: Path,
    config: dict[str, Any],
    run_id: str,
    device_str: str = "cuda",
    overfit_fixture: bool = False,
    study_holdout_fixture: bool = False,
    max_epochs: int | None = None,
    max_loci: int | None = None,
) -> TrainResult:
    """Train FlatDeepSet + tissue head; write run/checkpoint artifacts."""
    seed = int(config.get("experiment", {}).get("seed", 42))
    _set_seed(seed)

    train_cfg = config.get("training", {})
    fixture_mode = overfit_fixture or study_holdout_fixture
    if not fixture_mode and str(config.get("pilot", {}).get("mode", "")) == "study_holdout_fixture":
        study_holdout_fixture = True
        fixture_mode = True
    require_cuda = bool(train_cfg.get("require_cuda", False)) and not fixture_mode
    if fixture_mode and device_str.startswith("cuda") and not torch.cuda.is_available():
        device_str = "cpu"
        require_cuda = False
    device = resolve_device(device_str, require_cuda=require_cuda)
    use_amp = device.type == "cuda" and str(train_cfg.get("precision", "")).startswith("bf16")

    epochs = int(max_epochs if max_epochs is not None else train_cfg.get("max_epochs", 50))
    patience = int(train_cfg.get("early_stopping_patience", 10))
    lr = float(train_cfg.get("learning_rate", 1e-3))
    weight_decay = float(train_cfg.get("weight_decay", 1e-4))
    grad_clip = float(train_cfg.get("gradient_clip_norm", 2.0))
    batch_size = max(1, int(train_cfg.get("batch_size", 1)))
    model_cfg = config.get("model", {})

    fixture_records: list[FlatSampleRecord] | None = None
    train_phenotypes: list[SamplePhenotype] | None = None
    val_phenotypes: list[SamplePhenotype] | None = None
    test_phenotypes: list[SamplePhenotype] | None = None
    pilot_store: _PilotStore | None = None
    train_records: list[FlatSampleRecord] | None = None
    val_records: list[FlatSampleRecord] | None = None
    class_weights: torch.Tensor | None = None
    task_kind = "multiclass"
    age_mean = 0.0
    age_std = 1.0

    if study_holdout_fixture:
        task = str(config.get("pilot", {}).get("fixture_task", "tissue"))
        bundle = make_synthetic_study_holdout_bundle(seed=seed, task=task)
        records: list[FlatSampleRecord] = list(bundle["records"])
        class_names = list(bundle["class_names"])
        gene_ids = list(bundle["gene_ids"])
        n_genes = int(bundle["n_genes"])
        input_dim = int(bundle["input_dim"])
        n_classes = int(bundle["n_classes"])
        studies = list(bundle["studies"])
        split_manifest = build_study_grouped_split(
            bundle["sample_rows"],
            train_studies=studies[:-2],
            validation_studies=[studies[-2]],
            external_test_studies=[studies[-1]],
            split_id="study-holdout-fixture-v1",
        )
        train_ids = set(split_manifest["train_sample_ids"])
        val_ids = set(split_manifest["validation_sample_ids"])
        train_records = [r for r in records if r.sample_id in train_ids]
        val_records = [r for r in records if r.sample_id in val_ids]
        split = split_manifest
        class_weights = _class_weights([r.class_index for r in train_records], n_classes)
        if torch.allclose(class_weights, torch.ones_like(class_weights)):
            class_weights = None
    elif overfit_fixture:
        bundle = make_synthetic_overfit_bundle(seed=seed)
        fixture_records = list(bundle["records"])
        class_names = list(bundle["class_names"])
        gene_ids = list(bundle["gene_ids"])
        n_genes = int(bundle["n_genes"])
        input_dim = int(bundle["input_dim"])
        n_classes = int(bundle["n_classes"])
        train_records = fixture_records
        val_records = fixture_records
        split = {
            "mode": "overfit_fixture",
            "train_sample_ids": [r.sample_id for r in fixture_records],
            "val_sample_ids": [r.sample_id for r in fixture_records],
        }
        class_weights = _class_weights([r.class_index for r in fixture_records], n_classes)
        if torch.allclose(class_weights, torch.ones_like(class_weights)):
            class_weights = None
        task_kind = "multiclass"
        age_mean = 0.0
        age_std = 1.0
        test_phenotypes: list[SamplePhenotype] | None = None
    else:
        pilot = config.get("pilot", {})
        matrix_id = str(pilot["matrix_id"])
        graph_id = str(pilot["graph_id"])
        feature_set = str(
            pilot.get("static_feature_set") or config.get("features", {}).get("static_feature_set")
        )
        matrix_paths = matrix_store_paths(data_root / "canonical" / "matrices" / matrix_id)
        sample_index = read_sample_index(matrix_paths.sample_index_path)
        locus_index = read_locus_index(matrix_paths.locus_index_path)
        sample_ids = sample_index.sort_values("row_index")["sample_id"].astype(str).tolist()
        mode = str(pilot.get("mode", "gse35069"))
        task_kind = str(pilot.get("task", "multiclass"))
        test_phenotypes = None
        age_mean = 0.0
        age_std = 1.0

        if mode in {"multitask_hub", "deeprvat_hub"}:
            task_kind = "multitask"
            data_cfg = config.get("data", {})
            table_rel = Path(
                str(
                    pilot.get("sample_phenotype_table")
                    or data_cfg.get(
                        "sample_phenotype_table",
                        "canonical/phenotypes/sample_phenotype_table.parquet",
                    )
                )
            )
            table_path = table_rel if table_rel.is_absolute() else data_root / table_rel
            if not table_path.is_file():
                raise FileNotFoundError(
                    f"deeprvat/multitask hub requires sample_phenotype_table: {table_path}"
                )
            ont_rel = Path(
                str(
                    pilot.get("tissue_ontology")
                    or data_cfg.get(
                        "tissue_ontology",
                        "canonical/phenotypes/tissue_ontology.yaml",
                    )
                )
            )
            ont_path = ont_rel if ont_rel.is_absolute() else data_root / ont_rel
            if not ont_path.is_file():
                raise FileNotFoundError(
                    f"deeprvat/multitask hub requires tissue_ontology: {ont_path}"
                )
            ontology = load_tissue_ontology(ont_path)
            phenotypes, class_names = load_multitask_phenotypes(
                table_path,
                sample_ids=sample_ids,
                class_names=ontology.class_names,
            )
            sample_rows = [
                {
                    "sample_id": p.sample_id,
                    "study_id": p.study_id or p.donor_id,
                    "platform": p.platform,
                }
                for p in phenotypes
            ]
            auto_split = bool(pilot.get("auto_split", False))
            train_studies = [str(x) for x in pilot.get("train_studies", [])]
            if auto_split or not train_studies:
                split = partition_studies_by_sample_count(
                    sample_rows,
                    seed=seed,
                    train_fraction=float(config.get("splits", {}).get("train_fraction", 0.7)),
                    val_fraction=float(config.get("splits", {}).get("val_fraction", 0.15)),
                    split_id=str(pilot.get("split_id", f"{matrix_id}-auto-v1")),
                )
            else:
                split = build_study_grouped_split(
                    sample_rows,
                    train_studies=train_studies,
                    validation_studies=[str(x) for x in pilot["validation_studies"]],
                    external_test_studies=[str(x) for x in pilot.get("external_test_studies", [])],
                    split_id=str(pilot.get("split_id", f"{matrix_id}-multitask-v1")),
                )
            train_ids = set(split["train_sample_ids"])
            val_ids = set(split["validation_sample_ids"])
            test_ids = set(split.get("external_test_sample_ids") or [])
            train_phenotypes = [p for p in phenotypes if p.sample_id in train_ids]
            val_phenotypes = [p for p in phenotypes if p.sample_id in val_ids]
            test_phenotypes = [p for p in phenotypes if p.sample_id in test_ids]
            ages = [float(p.age) for p in train_phenotypes if p.age_mask and p.age is not None]
            if ages:
                age_mean = float(np.mean(ages))
                age_std = float(np.std(ages))
                if age_std < 1e-6:
                    age_std = 1.0
            else:
                age_mean = 0.0
                age_std = 1.0
        elif mode == "hub_pack":
            pheno_path = matrix_paths.root / "sample_phenotypes.parquet"
            if not pheno_path.is_file():
                raise FileNotFoundError(
                    f"hub_pack mode requires sample_phenotypes.parquet next to matrix: {pheno_path}"
                )
            value_column = str(pilot.get("label_column", "phenotype_value"))
            empty_as_control = bool(pilot.get("empty_as_control", False))
            if task_kind == "regression":
                phenotypes, class_names = load_hub_regression_phenotypes(
                    pheno_path, sample_ids=sample_ids
                )
            else:
                phenotypes, class_names = load_hub_sample_info_phenotypes(
                    pheno_path,
                    sample_ids=sample_ids,
                    value_column=value_column,
                    empty_as_control=empty_as_control,
                )
            sample_rows = [
                {
                    "sample_id": p.sample_id,
                    "study_id": p.study_id or p.donor_id,
                    "platform": p.platform,
                }
                for p in phenotypes
            ]
            split = build_study_grouped_split(
                sample_rows,
                train_studies=[str(x) for x in pilot["train_studies"]],
                validation_studies=[str(x) for x in pilot["validation_studies"]],
                external_test_studies=[str(x) for x in pilot.get("external_test_studies", [])],
                split_id=str(pilot.get("split_id", f"{matrix_id}-study-grouped-v1")),
            )
            train_ids = set(split["train_sample_ids"])
            val_ids = set(split["validation_sample_ids"])
            test_ids = set(split.get("external_test_sample_ids") or [])
            train_phenotypes = [p for p in phenotypes if p.sample_id in train_ids]
            val_phenotypes = [p for p in phenotypes if p.sample_id in val_ids]
            test_phenotypes = [p for p in phenotypes if p.sample_id in test_ids]
            if task_kind == "regression":
                ages = [float(p.age) for p in train_phenotypes if p.age is not None]
                if not ages:
                    raise ValueError("no train ages for standardization")
                age_mean = float(np.mean(ages))
                age_std = float(np.std(ages))
                if age_std < 1e-6:
                    age_std = 1.0
        else:
            meta_rel = Path(str(pilot["phenotype_metadata"]))
            metadata_path = meta_rel if meta_rel.is_absolute() else data_root / meta_rel
            phenotypes, class_names = load_gse35069_phenotypes(metadata_path, sample_ids=sample_ids)
            train_donors = {str(x) for x in pilot.get("train_donors", ["1", "2", "3", "4"])}
            val_donors = {str(x) for x in pilot.get("val_donors", ["5", "6"])}
            train_phenotypes, val_phenotypes, split = _split_by_donor(
                phenotypes, train_donors=train_donors, val_donors=val_donors
            )
            task_kind = "multiclass"

        lr_edges, regions = load_graph_tables(data_root / "canonical" / "graphs" / graph_id)
        locus_gene = build_locus_gene_index(
            locus_index=locus_index,
            locus_region_edges=lr_edges,
            regions=regions,
            max_loci=max_loci,
        )
        gene_ids = locus_gene.gene_ids
        n_genes = locus_gene.n_genes
        n_classes = 1 if task_kind == "regression" else len(class_names)
        n_cols = locus_gene.n_study_loci if max_loci is None else int(max_loci)

        static_paths = static_feature_store_paths(
            data_root / "canonical" / "static_features" / feature_set
        )
        static_by_col, static_valid, static_dim = build_static_column_table(
            locus_index_locus_ids=locus_index["locus_id"].to_numpy(),
            static_loci=read_loci_index(static_paths.loci_path),
            embeddings=open_embeddings_zarr(static_paths.embeddings_path),
            n_study_loci=n_cols,
        )
        epsilon = float(config.get("features", {}).get("methylation", {}).get("epsilon", 0.001))
        input_dim = 2 + static_dim
        sample_row_by_id = {
            str(sid): int(row)
            for sid, row in zip(
                sample_index["sample_id"].astype(str),
                sample_index["row_index"].astype(int),
                strict=True,
            )
        }
        # Prefer phenotype-table row_index when multitask (must match matrix).
        if mode in {"multitask_hub", "deeprvat_hub"}:
            table_rows = {str(p.sample_id): sample_row_by_id[str(p.sample_id)] for p in phenotypes}
            sample_row_by_id = table_rows
        pilot_store = _PilotStore(
            phenotypes=phenotypes,
            sample_row_by_id=sample_row_by_id,
            betas=open_betas_zarr(matrix_paths.betas_path),
            static_by_col=static_by_col,
            static_valid=static_valid,
            locus_gene=locus_gene,
            epsilon=epsilon,
            n_cols=n_cols,
        )
        if task_kind != "regression":
            tissue_train = [
                p.class_index for p in train_phenotypes if task_kind != "multitask" or p.tissue_mask
            ]
            if tissue_train:
                class_weights = _class_weights(tissue_train, n_classes)
                if torch.allclose(class_weights, torch.ones_like(class_weights)):
                    class_weights = None
            else:
                class_weights = None
        split.update(
            {
                "n_genes": n_genes,
                "n_classes": n_classes,
                "class_names": class_names,
                "matrix_id": matrix_id,
                "max_loci": max_loci,
                "task": task_kind,
                "age_mean": age_mean,
                "age_std": age_std,
            }
        )

    if study_holdout_fixture:
        task_kind = "multiclass"
        age_mean = 0.0
        age_std = 1.0
        test_phenotypes = None

    pool_name: PoolName = str(model_cfg.get("pooling", "max"))  # type: ignore[assignment]
    model = FlatDeepSet(
        input_dim,
        phi_hidden_dim=int(model_cfg.get("phi_hidden_dimension", 20)),
        phi_layers=int(model_cfg.get("phi_layers", 2)),
        rho_hidden_dim=int(model_cfg.get("rho_hidden_dimension", 10)),
        rho_layers=int(model_cfg.get("rho_layers", 3)),
        pool=pool_name,
        neutral_score=float(model_cfg.get("neutral_score", 0.5)),
        dropout=float(model_cfg.get("dropout", 0.0)),
    ).to(device)
    if task_kind == "regression":
        head: nn.Module = nn.Linear(n_genes, 1).to(device)
        with torch.no_grad():
            head.weight.normal_(0.0, 0.05)
            head.bias.zero_()
    elif task_kind == "multitask":
        _heads = config.get("heads", {})
        sex_cfg = _heads.get("sex", {}) if isinstance(_heads, dict) else {}
        sex_on = bool(sex_cfg.get("enabled", False))
        head = MultitaskHeads(n_genes, n_classes, sex_enabled=sex_on).to(device)
    else:
        seed_mask = torch.ones(n_classes, n_genes, dtype=torch.float32, device=device)
        head = SeedMaskedLinearHead(n_genes, n_classes, seed_mask).to(device)
        # Zero gene_weight is a saddle for CE on saturated MBS; break symmetry.
        with torch.no_grad():
            head.gene_weight.normal_(0.0, 0.05)

    loss_cfg = config.get("loss", {})
    lambda_age = float(loss_cfg.get("lambda_age", 1.0))
    lambda_tissue = float(loss_cfg.get("lambda_tissue", 1.0))
    lambda_sex = float(loss_cfg.get("lambda_sex", 1.0))
    heads_cfg = config.get("heads", {})
    age_head_cfg = heads_cfg.get("age", {}) if isinstance(heads_cfg, dict) else {}
    huber_delta = float(age_head_cfg.get("huber_delta", 1.0))
    age_loss_name = str(age_head_cfg.get("loss", "huber"))
    if age_loss_name not in {"huber", "mse"}:
        age_loss_name = "huber"
    opt_name = "adam" if overfit_fixture else str(train_cfg.get("optimizer", "adamw")).lower()
    if opt_name == "adam":
        optimizer = torch.optim.Adam(
            list(model.parameters()) + list(head.parameters()),
            lr=lr,
            weight_decay=weight_decay,
        )
    else:
        optimizer = torch.optim.AdamW(
            list(model.parameters()) + list(head.parameters()),
            lr=lr,
            weight_decay=weight_decay,
        )

    history: list[dict[str, Any]] = []
    best_val = float("inf")
    best_epoch = 0
    stale = 0
    cfg_hash = config_sha256(config)
    run_root = run_dir(artifact_root, run_id)
    ckpt_root = checkpoint_dir(artifact_root, run_id)
    run_root.mkdir(parents=True, exist_ok=True)
    ckpt_root.mkdir(parents=True, exist_ok=True)
    checkpoint_hashes: dict[str, str] = {}

    log_cfg = config.get("logging", {})
    use_tb = bool(log_cfg.get("tensorboard", False))
    # Default on when TensorBoard logging is enabled.
    auto_tb = bool(log_cfg.get("auto_tensorboard", use_tb))
    tb_port = int(log_cfg.get("tensorboard_port", 6006))
    tb_writer = None
    tb_server = None
    jsonl_path = run_root / "metrics.jsonl"
    if use_tb:
        from torch.utils.tensorboard import SummaryWriter  # noqa: PLC0415

        tb_dir = run_root / "tb"
        tb_dir.mkdir(parents=True, exist_ok=True)
        tb_writer = SummaryWriter(log_dir=str(tb_dir))
        if auto_tb:
            from mbs.training.monitor import ensure_tensorboard  # noqa: PLC0415

            try:
                tb_server = ensure_tensorboard(
                    run_root=run_root,
                    logdir=tb_dir,
                    preferred_port=tb_port,
                )
            except RuntimeError as tb_error:
                # Train must not die if TensorBoard is missing / port-stuck.
                warnings.warn(
                    f"auto TensorBoard skipped: {tb_error}",
                    stacklevel=2,
                )
                tb_server = None

    for epoch in range(1, epochs + 1):
        train_metrics = _run_epoch(
            records=train_records,
            phenotypes=train_phenotypes,
            pilot_store=pilot_store,
            model=model,
            head=head,
            optimizer=optimizer,
            device=device,
            n_genes=n_genes,
            class_weights=class_weights,
            use_amp=use_amp,
            grad_clip=grad_clip,
            train=True,
            task=task_kind,
            age_mean=age_mean,
            age_std=age_std,
            lambda_age=lambda_age,
            lambda_tissue=lambda_tissue,
            lambda_sex=lambda_sex,
            huber_delta=huber_delta,
            age_loss_name=age_loss_name,
            batch_size=batch_size,
        )
        val_metrics = _run_epoch(
            records=val_records,
            phenotypes=val_phenotypes,
            pilot_store=pilot_store,
            model=model,
            head=head,
            optimizer=None,
            device=device,
            n_genes=n_genes,
            class_weights=class_weights,
            use_amp=use_amp,
            grad_clip=grad_clip,
            train=False,
            task=task_kind,
            age_mean=age_mean,
            age_std=age_std,
            lambda_age=lambda_age,
            lambda_tissue=lambda_tissue,
            lambda_sex=lambda_sex,
            huber_delta=huber_delta,
            age_loss_name=age_loss_name,
            batch_size=batch_size,
        )
        row = {
            "epoch": epoch,
            "train_loss": train_metrics["loss"],
            "train_accuracy": train_metrics["accuracy"],
            "train_mae": train_metrics["mae"],
            "train_sex_accuracy": train_metrics.get("sex_accuracy", 0.0),
            "val_loss": val_metrics["loss"],
            "val_accuracy": val_metrics["accuracy"],
            "val_mae": val_metrics["mae"],
            "val_sex_accuracy": val_metrics.get("sex_accuracy", 0.0),
            "learning_rate": lr,
            "task": task_kind,
        }
        history.append(row)
        with jsonl_path.open("a", encoding="utf-8") as jsonl_handle:
            jsonl_handle.write(json.dumps(row) + "\n")
        if tb_writer is not None:
            tb_writer.add_scalar("loss/train", row["train_loss"], epoch)
            tb_writer.add_scalar("loss/val", row["val_loss"], epoch)
            if task_kind in {"regression", "multitask"}:
                tb_writer.add_scalar("mae/train", row["train_mae"], epoch)
                tb_writer.add_scalar("mae/val", row["val_mae"], epoch)
            if task_kind in {"multiclass", "multitask"}:
                tb_writer.add_scalar("accuracy/train", row["train_accuracy"], epoch)
                tb_writer.add_scalar("accuracy/val", row["val_accuracy"], epoch)
            if task_kind == "multitask":
                tb_writer.add_scalar("sex_accuracy/train", row["train_sex_accuracy"], epoch)
                tb_writer.add_scalar("sex_accuracy/val", row["val_sex_accuracy"], epoch)
            tb_writer.add_scalar("lr", lr, epoch)
        checkpoint_hashes["last.pt"] = save_checkpoint(
            ckpt_root / "last.pt",
            model_state=model.state_dict(),
            head_state=head.state_dict(),
            optimizer_state=optimizer.state_dict(),
            epoch=epoch,
            metrics=row,
            config_hash=cfg_hash,
        )

        if overfit_fixture:
            improved = train_metrics["accuracy"] >= 0.999 or train_metrics["loss"] < best_val - 1e-8
            best_val = min(best_val, train_metrics["loss"])
        else:
            improved = val_metrics["loss"] < best_val - 1e-6
            if improved:
                best_val = val_metrics["loss"]

        if improved:
            best_epoch = epoch
            stale = 0
            checkpoint_hashes["best.pt"] = save_checkpoint(
                ckpt_root / "best.pt",
                model_state=model.state_dict(),
                head_state=head.state_dict(),
                optimizer_state=optimizer.state_dict(),
                epoch=epoch,
                metrics=row,
                config_hash=cfg_hash,
            )
        else:
            stale += 1

        if overfit_fixture and train_metrics["accuracy"] >= 0.999:
            break
        if not overfit_fixture and stale >= patience:
            break

    holdout_metrics: dict[str, Any] | None = None
    if (
        test_phenotypes
        and pilot_store is not None
        and not overfit_fixture
        and not study_holdout_fixture
    ):
        # Reload best checkpoint for external_test scoring.
        best_path = ckpt_root / "best.pt"
        if best_path.is_file():
            payload = torch.load(best_path, map_location=device, weights_only=False)
            model.load_state_dict(payload["model_state"])
            head.load_state_dict(payload["head_state"])
        test_metrics = _run_epoch(
            records=None,
            phenotypes=test_phenotypes,
            pilot_store=pilot_store,
            model=model,
            head=head,
            optimizer=None,
            device=device,
            n_genes=n_genes,
            class_weights=class_weights,
            use_amp=use_amp,
            grad_clip=grad_clip,
            train=False,
            task=task_kind,
            age_mean=age_mean,
            age_std=age_std,
            lambda_age=lambda_age,
            lambda_tissue=lambda_tissue,
            lambda_sex=lambda_sex,
            huber_delta=huber_delta,
            age_loss_name=age_loss_name,
            batch_size=batch_size,
        )
        test_mae = test_metrics["mae"]
        if task_kind in {"regression", "multitask"}:
            test_mae = test_mae * age_std
        holdout_metrics = {
            "n_samples": int(test_metrics["n_samples"]),
            "loss": test_metrics["loss"],
            "accuracy": test_metrics["accuracy"],
            "mae": test_mae,
            "mae_note": (
                "years (destandardized)" if task_kind in {"regression", "multitask"} else None
            ),
            "age_n": test_metrics.get("age_n"),
            "tissue_n": test_metrics.get("tissue_n"),
            "sex_n": test_metrics.get("sex_n"),
            "sex_accuracy": test_metrics.get("sex_accuracy"),
        }

    age_std_meta = None
    if task_kind in {"regression", "multitask"}:
        age_std_meta = {"mean": age_mean, "std": age_std}
    metrics_out: dict[str, Any] = {
        "history": history,
        "best_epoch": best_epoch,
        "best_val_loss": best_val,
        "final": history[-1] if history else {},
        "n_genes": n_genes,
        "n_classes": n_classes,
        "class_names": class_names,
        "gene_panel_size": len(gene_ids),
        "overfit_fixture": overfit_fixture,
        "device": str(device),
        "model_public_name": "deepMAT",
        "task": task_kind,
        "external_test": holdout_metrics,
        "age_standardization": age_std_meta,
    }
    if overfit_fixture and history:
        metrics_out["overfit_train_accuracy"] = history[-1]["train_accuracy"]
        metrics_out["overfit_ok"] = bool(history[-1]["train_accuracy"] >= 0.999)

    if tb_writer is not None:
        tb_writer.flush()
        tb_writer.close()

    resolved = dict(config)
    resolved["runtime"] = {
        "run_id": run_id,
        "device": str(device),
        "input_dim": input_dim,
        "n_genes": n_genes,
        "n_classes": n_classes,
        "overfit_fixture": overfit_fixture,
        "max_loci": max_loci,
        "project_root": str(project_root),
    }
    write_run_artifacts(
        run_root=run_root,
        ckpt_root=ckpt_root,
        config=resolved,
        environment=collect_environment(device=device),
        metrics=metrics_out,
        split=split,
        checkpoint_hashes=checkpoint_hashes,
        config_hash=cfg_hash,
    )
    (run_root / "metrics.json").write_text(
        json.dumps(metrics_out, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return TrainResult(
        run_id=run_id,
        run_dir=run_root,
        checkpoint_dir=ckpt_root,
        metrics=metrics_out,
        best_epoch=best_epoch,
        tensorboard_url=None if tb_server is None else tb_server.url,
        tensorboard_port=None if tb_server is None else tb_server.port,
        monitor_hint=(f"uv run mbs monitor --run-id {run_id}" if use_tb else None),
    )
