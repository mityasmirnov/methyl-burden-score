"""HierarchicalDeepSet training on the 5d DeepRVAT multitask cohort (Milestone 6)."""

from __future__ import annotations

import json
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import torch

from mbs.evaluation.annotation_slices import (
    compare_hierarchical_vs_flat,
    index_annotation_summary,
)
from mbs.evaluation.metrics import multiclass_metrics, regression_metrics
from mbs.evaluation.splits import build_study_grouped_split
from mbs.matrix.store import (
    matrix_store_paths,
    read_locus_index,
    read_sample_index,
)
from mbs.matrix.virtual_hub_store import open_betas_for_matrix
from mbs.models import HierarchicalDeepSet
from mbs.scoring.orientation import (
    accumulate_signed_gene_mean_m,
    orient_run_scores,
    score_manifest,
)
from mbs.segment_ops import PoolName
from mbs.static_features.store import (
    open_embeddings_zarr,
    read_loci_index,
    static_feature_store_paths,
)
from mbs.training.encoder_config import resolve_encoder
from mbs.training.features import build_static_column_table, cpg_input_dim
from mbs.training.hier_dataset import (
    HierBatch,
    HierSampleRecord,
    build_hier_sample,
    make_synthetic_hier_overfit_bundle,
    pack_hier_records_to_batch,
)
from mbs.training.level1_norm import (
    Level1NormParams,
    fit_level1_from_betas,
    persist_level1,
    resolve_level1_config,
)
from mbs.training.cascade_assign import build_cascade_assignment, gene_linked_col_index
from mbs.training.locus_region_gene import (
    REGION_TYPE_TO_ID,
    RESIDUAL_PANEL_ID,
    LocusRegionGeneIndex,
    build_locus_region_gene_index,
    locus_region_gene_col_filter,
    region_systems_from_arm,
    region_type_vocab,
)
from mbs.training.loop import (
    TrainResult,
    load_experiment_config,
    maybe_constrained_split,
    resolve_device,
    split_sample_rows,
    task_key,
)
from mbs.training.multitask import MultitaskHeads, masked_multitask_loss
from mbs.training.phenotype_table import load_tissue_ontology
from mbs.training.phenotypes import SamplePhenotype, load_multitask_phenotypes
from mbs.training.run_artifacts import (
    checkpoint_dir,
    collect_environment,
    config_sha256,
    run_dir,
    save_checkpoint,
    write_run_artifacts,
)
from mbs.training.sampler import iter_epoch_batches

# Re-export for CLI convenience
__all__ = ["load_experiment_config", "train_hierarchical_baseline"]


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


def _label_flags_for_record(
    *,
    ph: SamplePhenotype | None,
    age_mean: float,
    age_std: float,
) -> tuple[float | None, bool, bool, bool, int]:
    age_enabled = bool(ph.age_mask) if ph is not None else False
    tissue_enabled = bool(ph.tissue_mask) if ph is not None else False
    sex_enabled = bool(ph.sex_mask) if ph is not None else False
    sex_cls = int(ph.sex_class_index) if ph is not None else 0
    age_value = None
    if age_enabled:
        if ph is None or ph.age is None:
            raise RuntimeError("missing age for masked sample")
        age_value = (float(ph.age) - age_mean) / age_std
    return age_value, age_enabled, tissue_enabled, sex_enabled, sex_cls


@dataclass(slots=True)
class _HierPilotStore:
    phenotypes: list[SamplePhenotype]
    sample_row_by_id: dict[str, int]
    betas: Any
    static_by_col: np.ndarray
    static_valid: np.ndarray
    locus_region: LocusRegionGeneIndex
    epsilon: float
    n_cols: int


def _packed_hier_mbs(
    model: HierarchicalDeepSet,
    batch: HierBatch,
    *,
    include_mapped: bool = True,
    include_residual: bool = True,
) -> tuple[torch.Tensor, torch.Tensor]:
    """Return panel MBS = [gene_mbs…, residual_slot] and presence masks."""
    batch_size = len(batch.sample_ids)
    n_regions = batch_size * max(int(batch.n_regions), 1)
    n_genes = batch_size * max(int(batch.n_genes), 1)
    residual_features = batch.residual_features if include_residual else batch.residual_features[:0]
    residual_sample_index = (
        batch.residual_sample_index if include_residual else batch.residual_sample_index[:0]
    )
    cpg_features = batch.cpg_features if include_mapped else batch.cpg_features[:0]
    cpg_to_region = batch.cpg_to_region if include_mapped else batch.cpg_to_region[:0]
    output = model(
        cpg_features=cpg_features,
        cpg_to_region=cpg_to_region,
        region_type=batch.region_type,
        region_to_gene=batch.region_to_gene,
        n_regions=n_regions if batch.n_regions > 0 else 0,
        n_gene_instances=n_genes if batch.n_genes > 0 else 0,
        residual_features=residual_features,
        residual_sample_index=residual_sample_index,
        n_samples=batch_size,
    )
    if batch.n_genes > 0:
        gene_mbs = output["mbs"].view(batch_size, batch.n_genes)
        gene_present = output["present"].view(batch_size, batch.n_genes)
    else:
        gene_mbs = torch.zeros(batch_size, 0, device=output["residual_mbs"].device)
        gene_present = torch.zeros(batch_size, 0, dtype=torch.bool, device=gene_mbs.device)
    if not include_mapped and batch.n_genes > 0:
        gene_mbs = torch.full_like(gene_mbs, float(model.neutral_score))
        gene_present = torch.zeros_like(gene_present)
    residual_mbs = output["residual_mbs"].view(batch_size, 1)
    residual_present = output["residual_present"].view(batch_size, 1)
    if not include_residual:
        residual_mbs = torch.full_like(residual_mbs, float(model.neutral_score))
        residual_present = torch.zeros_like(residual_present)
    mbs = torch.cat([gene_mbs, residual_mbs], dim=1)
    present = torch.cat([gene_present, residual_present], dim=1)
    return mbs, present


def _run_hier_epoch(
    *,
    records: list[HierSampleRecord] | None,
    phenotypes: list[SamplePhenotype] | None,
    pilot_store: _HierPilotStore | None,
    locus_region: LocusRegionGeneIndex,
    model: HierarchicalDeepSet,
    head: MultitaskHeads,
    optimizer: torch.optim.Optimizer | None,
    device: torch.device,
    class_weights: torch.Tensor | None,
    use_amp: bool,
    grad_clip: float,
    train: bool,
    age_mean: float,
    age_std: float,
    lambda_age: float,
    lambda_tissue: float,
    lambda_sex: float,
    huber_delta: float,
    age_loss_name: str,
    batch_size: int,
    seed: int = 42,
    epoch: int = 0,
    batch_token_budget: int | None = None,
    allowed_region_type_ids: set[int] | None = None,
    include_residual: bool = True,
    include_mapped: bool = True,
    include_m_value: bool = True,
    include_robust_z: bool = False,
    level1_params: Level1NormParams | None = None,
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
    pred_tissue: list[int] = []
    true_tissue: list[int] = []
    pred_age: list[float] = []
    true_age: list[float] = []
    amp_dtype = torch.bfloat16 if device.type == "cuda" else torch.float32
    context = torch.enable_grad() if train else torch.no_grad()
    ph_by_id = {p.sample_id: p for p in (phenotypes or [])}

    if records is not None:
        all_records: list[HierSampleRecord] | None = records
        all_phenotypes: list[SamplePhenotype] | None = None
    else:
        all_records = None
        if phenotypes is None or pilot_store is None:
            raise RuntimeError("pilot phenotypes and store are required when records is None")
        all_phenotypes = phenotypes
    n_items = len(all_records) if all_records is not None else len(all_phenotypes or [])
    n_tokens: list[int] = []
    study_ids: list[str] = []
    task_keys: list[str] = []
    if all_records is not None:
        for rec in all_records:
            ph = ph_by_id.get(rec.sample_id)
            n_tokens.append(
                max(1, rec.features.n_observed_edges + rec.features.n_observed_residual)
            )
            study_ids.append(str(ph.study_id if ph and ph.study_id else rec.sample_id))
            task_keys.append(task_key(ph))
    else:
        for ph in all_phenotypes or []:
            n_tokens.append(1)
            study_ids.append(str(ph.study_id or ph.sample_id))
            task_keys.append(task_key(ph))

    with context:
        first_batch = True
        for idxs in iter_epoch_batches(
            n_items,
            n_tokens=n_tokens,
            study_ids=study_ids,
            task_keys=task_keys,
            batch_token_budget=batch_token_budget,
            batch_size=batch_size,
            seed=seed,
            epoch=epoch,
        ):
            if all_records is not None:
                chunk = [all_records[i] for i in idxs]
            else:
                if all_phenotypes is None or pilot_store is None:
                    raise RuntimeError("missing phenotypes/store for materialization")
                chunk = [
                    _materialize_hier_record(
                        all_phenotypes[i],
                        pilot_store,
                        include_m_value=include_m_value,
                        include_robust_z=include_robust_z,
                        level1_params=level1_params,
                    )
                    for i in idxs
                ]
            age_values: list[float | None] = []
            age_flags: list[bool] = []
            tissue_flags: list[bool] = []
            sex_flags: list[bool] = []
            sex_idxs: list[int] = []
            for record in chunk:
                ph = ph_by_id.get(record.sample_id)
                age_value, age_on, tissue_on, sex_on, sex_cls = _label_flags_for_record(
                    ph=ph,
                    age_mean=age_mean,
                    age_std=age_std,
                )
                age_values.append(age_value)
                age_flags.append(age_on)
                tissue_flags.append(tissue_on)
                sex_flags.append(sex_on)
                sex_idxs.append(sex_cls)

            batch = pack_hier_records_to_batch(
                chunk,
                locus_region=locus_region,
                age_values=age_values,
                age_enabled=age_flags,
                tissue_enabled=tissue_flags,
                sex_enabled=sex_flags,
                sex_class_indices=sex_idxs,
                allowed_region_type_ids=allowed_region_type_ids,
                include_residual=include_residual,
            ).to(device)

            if train and optimizer is not None:
                optimizer.zero_grad(set_to_none=True)
            with torch.autocast(device_type=device.type, dtype=amp_dtype, enabled=use_amp):
                mbs, present = _packed_hier_mbs(
                    model,
                    batch,
                    include_mapped=include_mapped,
                    include_residual=include_residual,
                )
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
            if first_batch:
                first_batch = False
                print(  # noqa: T201
                    f"[hier] first batch packed "
                    f"cpg={tuple(batch.cpg_features.shape)} "
                    f"residual={tuple(batch.residual_features.shape)} "
                    f"regions={batch.n_regions} genes={batch.n_genes}",
                    flush=True,
                )
            if train and optimizer is not None:
                loss.backward()
                if grad_clip > 0:
                    torch.nn.utils.clip_grad_norm_(
                        list(model.parameters()) + list(head.parameters()),
                        grad_clip,
                    )
                optimizer.step()

            total_loss += float(loss.detach().item()) * len(chunk)
            total_correct += float(result.metrics.get("tissue_correct", 0.0))
            tissue_n += float(result.metrics.get("tissue_n", 0.0))
            step_age_n = float(result.metrics.get("age_n", 0.0))
            total_mae_sum += float(result.metrics.get("mae", 0.0)) * step_age_n
            age_n += step_age_n
            sex_n += float(result.metrics.get("sex_n", 0.0))
            sex_correct += float(result.metrics.get("sex_correct", 0.0))
            n += len(chunk)
            if not train:
                # Autocast may leave activations in bf16; heads are fp32.
                mbs_f = mbs.float()
                present_f = present
                if bool(batch.tissue_mask.any()):
                    logits = head.forward_tissue(mbs_f, present_f)
                    tmask = batch.tissue_mask.reshape(-1)
                    pred_tissue.extend(
                        logits.argmax(dim=-1).reshape(-1)[tmask].detach().cpu().tolist()
                    )
                    true_tissue.extend(
                        batch.tissue_target.reshape(-1)[tmask].detach().cpu().tolist()
                    )
                if batch.age_target is not None and bool(batch.age_mask.any()):
                    age_hat = head.forward_age(mbs_f, present_f)
                    amask = batch.age_mask.reshape(-1)
                    pred_age.extend(age_hat.reshape(-1)[amask].detach().cpu().tolist())
                    true_age.extend(batch.age_target.reshape(-1)[amask].detach().cpu().tolist())

    out: dict[str, float] = {
        "loss": total_loss / max(n, 1),
        "accuracy": total_correct / max(tissue_n, 1.0),
        "mae": total_mae_sum / max(age_n, 1.0),
        "sex_accuracy": sex_correct / max(sex_n, 1.0),
        "n_samples": float(n),
        "age_n": age_n,
        "tissue_n": tissue_n,
        "sex_n": sex_n,
    }
    if pred_tissue and true_tissue:
        tissue_m = multiclass_metrics(np.asarray(true_tissue), np.asarray(pred_tissue))
        out["macro_f1"] = float(tissue_m.get("macro_f1") or 0.0)
        out["balanced_accuracy"] = float(tissue_m.get("balanced_accuracy") or 0.0)
    if pred_age and true_age:
        # Match flat loop: RMSE/R² on standardized age; MAE destandardized in holdout.
        age_m = regression_metrics(
            np.asarray(true_age, dtype=np.float64),
            np.asarray(pred_age, dtype=np.float64),
        )
        out["rmse"] = float(age_m.get("rmse") or 0.0)
        out["r2"] = float(age_m.get("r2") or 0.0)
        out["pearson_r"] = float(age_m.get("pearson_r") or 0.0)
        out["spearman_r"] = float(age_m.get("spearman_r") or 0.0)
    return out


def _materialize_hier_record(
    phenotype: SamplePhenotype,
    store: _HierPilotStore,
    *,
    include_m_value: bool = True,
    include_robust_z: bool = False,
    level1_params: Any | None = None,
) -> HierSampleRecord:
    row = store.sample_row_by_id[phenotype.sample_id]
    beta_row = np.asarray(store.betas[row, : store.n_cols], dtype=np.float32)
    return build_hier_sample(
        phenotype=phenotype,
        beta_row=beta_row,
        static_by_col=store.static_by_col,
        static_valid=store.static_valid,
        locus_region=store.locus_region,
        epsilon=store.epsilon,
        include_m_value=include_m_value,
        include_robust_z=include_robust_z,
        level1_params=level1_params,
    )


def _load_reuse_flat_split(
    artifact_root: Path,
    *,
    flat_run_id: str,
) -> dict[str, Any] | None:
    path = artifact_root / "runs" / flat_run_id / "split.json"
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def _promoter_type_ids() -> set[int]:
    return {
        REGION_TYPE_TO_ID["promoter_core"],
        REGION_TYPE_TO_ID["promoter_proximal"],
    }


def train_hierarchical_baseline(
    *,
    project_root: Path,
    data_root: Path,
    artifact_root: Path,
    config: dict[str, Any],
    run_id: str,
    device_str: str = "cuda",
    overfit_fixture: bool = False,
    max_epochs: int | None = None,
    max_loci: int | None = None,
) -> TrainResult:
    """Train HierarchicalDeepSet + MultitaskHeads on the 5d cohort (or fixture)."""
    seed = int(config.get("experiment", {}).get("seed", 42))
    _set_seed(seed)

    train_cfg = config.get("training", {})
    level1_cfg = resolve_level1_config(config)
    include_m_value = bool(level1_cfg["include_m_value"])
    include_robust_z = bool(level1_cfg["include_robust_z"])
    level1_epsilon = float(level1_cfg["epsilon"])
    level1_sigma_min = float(level1_cfg["sigma_min"])
    level1_params: Level1NormParams | None = None
    level1_manifest: dict[str, Any] | None = None
    require_cuda = bool(train_cfg.get("require_cuda", False)) and not overfit_fixture
    if overfit_fixture and device_str.startswith("cuda") and not torch.cuda.is_available():
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
    raw_budget = train_cfg.get("batch_token_budget")
    batch_token_budget = int(raw_budget) if raw_budget not in (None, "", 0) else None
    model_cfg = config.get("model", {})
    arm = model_cfg.get("arm")
    if model_cfg.get("region_systems"):
        region_systems = tuple(str(s) for s in model_cfg["region_systems"])
    else:
        region_systems = region_systems_from_arm(str(arm) if arm is not None else None)
    type_vocab = region_type_vocab(region_systems)

    train_records: list[HierSampleRecord] | None = None
    val_records: list[HierSampleRecord] | None = None
    train_phenotypes: list[SamplePhenotype] | None = None
    val_phenotypes: list[SamplePhenotype] | None = None
    test_phenotypes: list[SamplePhenotype] | None = None
    pilot_store: _HierPilotStore | None = None
    class_weights: torch.Tensor | None = None
    age_mean = 0.0
    age_std = 1.0

    if overfit_fixture:
        if region_systems != ("gene",):
            raise ValueError(
                "hier overfit_fixture supports gene-only; use flat branch arms for rbs/tbs fixtures"
            )
        bundle = make_synthetic_hier_overfit_bundle(
            seed=seed,
            include_m_value=include_m_value,
            include_robust_z=include_robust_z,
            sigma_min=level1_sigma_min,
            epsilon=level1_epsilon,
        )
        train_records = list(bundle["records"])
        val_records = train_records
        locus_region = bundle["locus_region"]
        class_names = list(bundle["class_names"])
        gene_ids = list(bundle["gene_ids"])
        n_genes = int(bundle["n_genes"])
        n_panel = int(bundle["n_panel"])
        input_dim = int(bundle["input_dim"])
        n_classes = int(bundle["n_classes"])
        level1_params = bundle.get("level1_params")
        split = {
            "mode": "overfit_fixture",
            "train_sample_ids": [r.sample_id for r in train_records],
            "val_sample_ids": [r.sample_id for r in train_records],
        }
        # Attach synthetic ages for multitask loss via phenotype stubs
        train_phenotypes = [
            SamplePhenotype(
                sample_id=r.sample_id,
                cell_type=class_names[r.class_index],
                donor_id=r.donor_id,
                title=r.sample_id,
                class_index=r.class_index,
                study_id="SYN",
                age=float(bundle["ages"][i]),
                platform="HM450",
                age_mask=True,
                tissue_mask=True,
                sex_mask=True,
                sex_class_index=i % 2,
            )
            for i, r in enumerate(train_records)
        ]
        val_phenotypes = train_phenotypes
        class_weights = _class_weights([r.class_index for r in train_records], n_classes)
        if torch.allclose(class_weights, torch.ones_like(class_weights)):
            class_weights = None
        ages = [float(a) for a in bundle["ages"]]
        age_mean = float(np.mean(ages))
        age_std = float(np.std(ages)) or 1.0
    else:
        pilot = config.get("pilot", {})
        mode = str(pilot.get("mode", "deeprvat_hub"))
        if mode not in {"multitask_hub", "deeprvat_hub"}:
            raise ValueError(
                "hierarchical train expects pilot.mode multitask_hub or deeprvat_hub "
                f"(got {mode!r})"
            )
        matrix_id = str(pilot["matrix_id"])
        graph_id = str(pilot["graph_id"])
        feature_set = str(
            pilot.get("static_feature_set") or config.get("features", {}).get("static_feature_set")
        )
        matrix_paths = matrix_store_paths(data_root / "canonical" / "matrices" / matrix_id)
        sample_index = read_sample_index(matrix_paths.sample_index_path)
        locus_index = read_locus_index(matrix_paths.locus_index_path)
        sample_ids = sample_index.sort_values("row_index")["sample_id"].astype(str).tolist()

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
        ont_rel = Path(
            str(
                pilot.get("tissue_ontology")
                or data_cfg.get("tissue_ontology", "canonical/phenotypes/tissue_ontology.yaml")
            )
        )
        ont_path = ont_rel if ont_rel.is_absolute() else data_root / ont_rel
        ontology = load_tissue_ontology(ont_path)
        phenotypes, class_names = load_multitask_phenotypes(
            table_path,
            sample_ids=sample_ids,
            class_names=ontology.class_names,
        )
        sample_rows = split_sample_rows(phenotypes)

        reuse = bool(pilot.get("reuse_flat_split", True))
        flat_run = str(
            pilot.get("flat_split_run_id", "stage0-flat-deeprvat-age-tissue-sex-full-v1")
        )
        split = _load_reuse_flat_split(artifact_root, flat_run_id=flat_run) if reuse else None
        if split is None:
            auto_split = bool(pilot.get("auto_split", False))
            train_studies = [str(x) for x in pilot.get("train_studies", [])]
            if auto_split or not train_studies:
                split = maybe_constrained_split(
                    sample_rows,
                    seed=seed,
                    train_fraction=float(config.get("splits", {}).get("train_fraction", 0.7)),
                    val_fraction=float(config.get("splits", {}).get("val_fraction", 0.15)),
                    split_id=str(pilot.get("split_id", f"{matrix_id}-hier-auto-v1")),
                )
            else:
                split = build_study_grouped_split(
                    sample_rows,
                    train_studies=train_studies,
                    validation_studies=[str(x) for x in pilot["validation_studies"]],
                    external_test_studies=[str(x) for x in pilot.get("external_test_studies", [])],
                    split_id=str(pilot.get("split_id", f"{matrix_id}-hier-v1")),
                )
            split["reused_flat_split"] = False
        else:
            split = dict(split)
            split["reused_flat_split"] = True
            split["flat_split_run_id"] = flat_run

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

        lr_edges, regions = load_graph_tables(data_root / "canonical" / "graphs" / graph_id)
        print(  # noqa: T201
            f"[hier] building locus→region index (max_loci={max_loci})",
            flush=True,
        )
        locus_region = build_locus_region_gene_index(
            locus_index=locus_index,
            locus_region_edges=lr_edges,
            regions=regions,
            max_loci=max_loci,
            region_systems=region_systems,
        )
        gene_linked_only = bool(train_cfg.get("gene_linked_only", False))
        if gene_linked_only:
            genes_path = data_root / "canonical" / "graphs" / graph_id / "genes.parquet"
            genes_df = (
                pd.read_parquet(genes_path) if genes_path.is_file() else pd.DataFrame()
            )
            cascade_assignment = build_cascade_assignment(
                locus_index=locus_index,
                locus_region_edges=lr_edges,
                regions=regions,
                genes=genes_df,
                max_loci=max_loci,
            )
            gene_cols = gene_linked_col_index(cascade_assignment)
            locus_region = locus_region_gene_col_filter(locus_region, gene_cols)
            print(  # noqa: T201
                f"[hier] gene_linked_only panel: {gene_cols.size} CpG columns, "
                f"{locus_region.n_typed_edges} typed edges",
                flush=True,
            )
        gene_ids = locus_region.gene_ids
        n_genes = locus_region.n_genes
        n_panel = locus_region.n_panel
        n_classes = len(class_names)
        n_cols = locus_region.n_study_loci if max_loci is None else int(max_loci)
        print(  # noqa: T201
            "[hier] index ready: "
            f"genes={n_genes} regions={locus_region.n_regions} "
            f"typed_edges={locus_region.n_typed_edges} "
            f"residual_cols={locus_region.n_residual_cols} study_loci={n_cols}",
            flush=True,
        )

        use_cpgpt = bool(config.get("stage0", {}).get("use_cpgpt_static_features", True))
        if not use_cpgpt or not feature_set or str(feature_set).lower() in {"none", "null", "off"}:
            static_by_col = np.zeros((n_cols, 0), dtype=np.float32)
            static_valid = np.zeros(n_cols, dtype=bool)
            static_dim = 0
            print("[hier] CpGPT static disabled (dim=0)", flush=True)  # noqa: T201
        else:
            static_paths = static_feature_store_paths(
                data_root / "canonical" / "static_features" / feature_set
            )
            print("[hier] aligning static features…", flush=True)  # noqa: T201
            static_by_col, static_valid, static_dim = build_static_column_table(
                locus_index_locus_ids=locus_index["locus_id"].to_numpy(),
                static_loci=read_loci_index(static_paths.loci_path),
                embeddings=open_embeddings_zarr(static_paths.embeddings_path),
                n_study_loci=n_cols,
            )
            print(  # noqa: T201
                f"[hier] static ready dim={static_dim} "
                f"valid_cols={int(static_valid.sum())}/{n_cols}",
                flush=True,
            )
        epsilon = float(level1_epsilon)
        input_dim = cpg_input_dim(
            static_dim,
            include_m_value=include_m_value,
            include_robust_z=include_robust_z,
        )
        sample_row_by_id = {
            str(sid): int(row)
            for sid, row in zip(
                sample_index["sample_id"].astype(str),
                sample_index["row_index"].astype(int),
                strict=True,
            )
        }
        sample_row_by_id = {
            str(p.sample_id): sample_row_by_id[str(p.sample_id)] for p in phenotypes
        }
        pilot_store = _HierPilotStore(
            phenotypes=phenotypes,
            sample_row_by_id=sample_row_by_id,
            betas=open_betas_for_matrix(matrix_paths.root),
            static_by_col=static_by_col,
            static_valid=static_valid,
            locus_region=locus_region,
            epsilon=epsilon,
            n_cols=n_cols,
        )
        if include_robust_z:
            if not train_phenotypes:
                raise ValueError("robust_deviation requires train phenotypes for Hub/pilot path")
            train_rows = [sample_row_by_id[str(p.sample_id)] for p in train_phenotypes]
            level1_params = fit_level1_from_betas(
                pilot_store.betas,
                train_rows,
                epsilon=level1_epsilon,
                sigma_min=level1_sigma_min,
                n_loci=n_cols,
                locus_ids=locus_index["locus_id"].to_numpy()[:n_cols],
                fold_id=str(split.get("split_id") or run_id),
                run_id=run_id,
            )
        tissue_train = [p.class_index for p in train_phenotypes if p.tissue_mask]
        if tissue_train:
            class_weights = _class_weights(tissue_train, n_classes)
            if torch.allclose(class_weights, torch.ones_like(class_weights)):
                class_weights = None
        split.update(
            {
                "n_genes": n_genes,
                "n_panel": n_panel,
                "n_regions": locus_region.n_regions,
                "n_classes": n_classes,
                "class_names": class_names,
                "matrix_id": matrix_id,
                "max_loci": max_loci,
                "task": "multitask",
                "age_mean": age_mean,
                "age_std": age_std,
                "n_typed_edges": locus_region.n_typed_edges,
                "n_residual_cols": locus_region.n_residual_cols,
                "annotation_summary": index_annotation_summary(locus_region),
                "region_types": list(locus_region.region_types),
                "region_systems": list(locus_region.region_systems),
            }
        )

    region_types_cfg = model_cfg.get("region_types") or list(type_vocab)
    if list(region_types_cfg) != list(type_vocab):
        raise ValueError(
            f"model.region_types must match {list(type_vocab)}, got {region_types_cfg}"
        )
    n_region_types = len(type_vocab)
    # Prefer index vocab when real graph was loaded
    if not overfit_fixture:
        n_region_types = len(locus_region.region_types)
        type_vocab = locus_region.region_types
    enc = resolve_encoder(
        model_cfg,
        default_activation="gelu",
        default_layer_norm=True,
        default_dropout=0.1,
        default_cpg_hidden=64,
    )
    cpg_pool: PoolName = str(model_cfg.get("cpg_pooling", "max"))  # type: ignore[assignment]
    region_pool: PoolName = str(model_cfg.get("region_pooling", "max"))  # type: ignore[assignment]
    model = HierarchicalDeepSet(
        input_dim,
        n_region_types,
        cpg_hidden_dim=int(enc["cpg_hidden_dim"]),
        region_hidden_dim=int(model_cfg.get("region_hidden_dimension", 32)),
        region_type_dim=int(model_cfg.get("region_type_dimension", 8)),
        cpg_pool=cpg_pool,
        region_pool=region_pool,
        residual_pool=str(model_cfg.get("residual_pooling", "max")),  # type: ignore[arg-type]
        neutral_score=float(model_cfg.get("neutral_score", 0.5)),
        dropout=float(enc["dropout"]),
        activation=str(enc["activation"]),
        layer_norm=bool(enc["layer_norm"]),
    ).to(device)

    heads_cfg = config.get("heads", {})
    sex_cfg = heads_cfg.get("sex", {}) if isinstance(heads_cfg, dict) else {}
    sex_on = bool(sex_cfg.get("enabled", True))
    head = MultitaskHeads(n_panel, n_classes, sex_enabled=sex_on).to(device)

    loss_cfg = config.get("loss", {})
    lambda_age = float(loss_cfg.get("lambda_age", 1.0))
    lambda_tissue = float(loss_cfg.get("lambda_tissue", 1.0))
    lambda_sex = float(loss_cfg.get("lambda_sex", 1.0))
    age_head_cfg = heads_cfg.get("age", {}) if isinstance(heads_cfg, dict) else {}
    huber_delta = float(age_head_cfg.get("huber_delta", 1.0))
    age_loss_name = str(age_head_cfg.get("loss", "huber"))
    if age_loss_name not in {"huber", "mse"}:
        age_loss_name = "huber"

    opt_name = "adam" if overfit_fixture else str(train_cfg.get("optimizer", "adamw")).lower()
    params = list(model.parameters()) + list(head.parameters())
    if opt_name == "adam":
        optimizer = torch.optim.Adam(params, lr=lr, weight_decay=weight_decay)
    else:
        optimizer = torch.optim.AdamW(params, lr=lr, weight_decay=weight_decay)

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
    if level1_params is not None and include_robust_z:
        if level1_params.run_id is None:
            level1_params = Level1NormParams(
                mu=level1_params.mu,
                sigma=level1_params.sigma,
                estimated=level1_params.estimated,
                locus_ids=level1_params.locus_ids,
                sigma_min=level1_params.sigma_min,
                n_train_samples=level1_params.n_train_samples,
                epsilon=level1_params.epsilon,
                fold_id=level1_params.fold_id or run_id,
                run_id=run_id,
            )
        level1_manifest = persist_level1(run_root, level1_params)

    log_cfg = config.get("logging", {})
    use_tb = bool(log_cfg.get("tensorboard", False))
    auto_tb = bool(log_cfg.get("auto_tensorboard", use_tb))
    tb_port = int(log_cfg.get("tensorboard_port", 6008))
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
                warnings.warn(f"auto TensorBoard skipped: {tb_error}", stacklevel=2)
                tb_server = None

    def run_epoch(
        *,
        records: list[HierSampleRecord] | None,
        phenotypes: list[SamplePhenotype] | None,
        optimizer: torch.optim.Optimizer | None,
        train: bool,
        epoch: int = 0,
        allowed_region_type_ids: set[int] | None = None,
        include_residual: bool = True,
        include_mapped: bool = True,
    ) -> dict[str, float]:
        return _run_hier_epoch(
            records=records,
            phenotypes=phenotypes,
            pilot_store=pilot_store,
            locus_region=locus_region,
            model=model,
            head=head,
            optimizer=optimizer,
            device=device,
            class_weights=class_weights,
            use_amp=use_amp,
            grad_clip=grad_clip,
            train=train,
            age_mean=age_mean,
            age_std=age_std,
            lambda_age=lambda_age,
            lambda_tissue=lambda_tissue,
            lambda_sex=lambda_sex,
            huber_delta=huber_delta,
            age_loss_name=age_loss_name,
            batch_size=batch_size,
            seed=seed,
            epoch=epoch,
            batch_token_budget=batch_token_budget,
            allowed_region_type_ids=allowed_region_type_ids,
            include_residual=include_residual,
            include_mapped=include_mapped,
            include_m_value=include_m_value,
            include_robust_z=include_robust_z,
            level1_params=level1_params,
        )

    print(  # noqa: T201
        f"[hier] training start epochs={epochs} batch_size={batch_size} "
        f"n_train={len(train_phenotypes or train_records or [])} "
        f"n_val={len(val_phenotypes or val_records or [])} device={device}",
        flush=True,
    )
    for epoch in range(1, epochs + 1):
        print(f"[hier] epoch {epoch}/{epochs} train…", flush=True)  # noqa: T201
        train_metrics = run_epoch(
            records=train_records,
            phenotypes=train_phenotypes,
            optimizer=optimizer,
            train=True,
            epoch=epoch,
        )
        val_metrics = run_epoch(
            records=val_records,
            phenotypes=val_phenotypes,
            optimizer=None,
            train=False,
            epoch=epoch,
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
            "task": "multitask",
        }
        history.append(row)
        with jsonl_path.open("a", encoding="utf-8") as jsonl_handle:
            jsonl_handle.write(json.dumps(row) + "\n")
        if tb_writer is not None:
            tb_writer.add_scalar("loss/train", row["train_loss"], epoch)
            tb_writer.add_scalar("loss/val", row["val_loss"], epoch)
            tb_writer.add_scalar("mae/train", row["train_mae"], epoch)
            tb_writer.add_scalar("mae/val", row["val_mae"], epoch)
            tb_writer.add_scalar("accuracy/train", row["train_accuracy"], epoch)
            tb_writer.add_scalar("accuracy/val", row["val_accuracy"], epoch)
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

    # Reload best for holdout + ablations
    best_path = ckpt_root / "best.pt"
    if best_path.is_file():
        payload = torch.load(best_path, map_location=device, weights_only=False)
        model.load_state_dict(payload["model_state"])
        head.load_state_dict(payload["head_state"])

    holdout_metrics: dict[str, Any] | None = None
    ablations: dict[str, Any] = {}
    eval_phenotypes = test_phenotypes if test_phenotypes else val_phenotypes
    eval_records = None if pilot_store is not None else val_records
    ablation_cap = int(config.get("evaluation", {}).get("ablation_max_samples", 512))
    if ablation_cap > 0 and eval_phenotypes:
        ablation_phenotypes = eval_phenotypes[:ablation_cap]
    else:
        ablation_phenotypes = eval_phenotypes
    if ablation_cap > 0 and eval_records:
        ablation_records = eval_records[:ablation_cap]
    else:
        ablation_records = eval_records
    if eval_phenotypes or eval_records:
        base = run_epoch(
            records=eval_records,
            phenotypes=eval_phenotypes if pilot_store is not None else train_phenotypes,
            optimizer=None,
            train=False,
        )
        mae = base["mae"] * age_std
        holdout_metrics = {
            "n_samples": int(base["n_samples"]),
            "loss": base["loss"],
            "accuracy": base["accuracy"],
            "mae": mae,
            "mae_note": "years (destandardized)",
            "age_n": base.get("age_n"),
            "tissue_n": base.get("tissue_n"),
            "sex_n": base.get("sex_n"),
            "sex_accuracy": base.get("sex_accuracy"),
            "macro_f1": base.get("macro_f1"),
            "balanced_accuracy": base.get("balanced_accuracy"),
            "rmse": base.get("rmse"),
            "r2": base.get("r2"),
            "pearson_r": base.get("pearson_r"),
            "spearman_r": base.get("spearman_r"),
            "split": "external_test" if test_phenotypes else "validation",
        }
        body_id = REGION_TYPE_TO_ID["gene_body"]
        ablation_specs: dict[str, dict[str, Any]] = {
            "full": {},
            "mapped_only": {"include_residual": False},
            "residual_only": {"include_mapped": False, "include_residual": True},
            "promoters_only": {
                "allowed_region_type_ids": _promoter_type_ids(),
                "include_residual": False,
            },
            "gene_body_only": {
                "allowed_region_type_ids": {body_id},
                "include_residual": False,
            },
        }
        for name, extra in ablation_specs.items():
            try:
                m = run_epoch(
                    records=ablation_records if pilot_store is None else None,
                    phenotypes=(
                        ablation_phenotypes if pilot_store is not None else train_phenotypes
                    ),
                    optimizer=None,
                    train=False,
                    allowed_region_type_ids=extra.get("allowed_region_type_ids"),
                    include_residual=bool(extra.get("include_residual", True)),
                    include_mapped=bool(extra.get("include_mapped", True)),
                )
            except ValueError as ablation_error:
                ablations[name] = {"skipped": True, "error": str(ablation_error)}
                continue
            ablations[name] = {
                "loss": m["loss"],
                "accuracy": m["accuracy"],
                "mae": m["mae"] * age_std,
                "sex_accuracy": m.get("sex_accuracy"),
                "n_samples": int(m["n_samples"]),
                "ablation_max_samples": ablation_cap,
                "slice": name,
            }

    flat_metrics = None
    flat_run_id = str(
        config.get("pilot", {}).get(
            "flat_split_run_id",
            "stage0-flat-deeprvat-age-tissue-sex-full-v1",
        )
    )
    flat_metrics_path = artifact_root / "runs" / flat_run_id / "metrics.json"
    if flat_metrics_path.is_file():
        try:
            flat_payload = json.loads(flat_metrics_path.read_text(encoding="utf-8"))
            flat_metrics = flat_payload.get("external_test") or flat_payload.get("final")
        except (OSError, json.JSONDecodeError):
            flat_metrics = None

    metrics_out: dict[str, Any] = {
        "history": history,
        "best_epoch": best_epoch,
        "best_val_loss": best_val,
        "final": history[-1] if history else {},
        "n_genes": n_genes,
        "n_panel": n_panel,
        "n_regions": locus_region.n_regions,
        "n_typed_edges": locus_region.n_typed_edges,
        "n_residual_cols": locus_region.n_residual_cols,
        "annotation_summary": index_annotation_summary(locus_region),
        "n_classes": n_classes,
        "class_names": class_names,
        "gene_panel_size": len(gene_ids),
        "panel_ids": [*gene_ids, RESIDUAL_PANEL_ID],
        "region_types": list(locus_region.region_types),
        "region_systems": list(locus_region.region_systems),
        "overfit_fixture": overfit_fixture,
        "device": str(device),
        "model_public_name": "deepMAT-hierarchical",
        "task": "multitask",
        "external_test": holdout_metrics,
        "ablations": ablations,
        "annotation_slices": {
            name: ablations[name]
            for name in ("full", "mapped_only", "residual_only")
            if name in ablations
        },
        "vs_flat": compare_hierarchical_vs_flat(
            hierarchical_metrics=holdout_metrics or {},
            flat_metrics=flat_metrics if isinstance(flat_metrics, dict) else None,
        ),
        "age_standardization": {"mean": age_mean, "std": age_std},
        "level1_normalization": (
            {
                "enabled": True,
                "manifest_hash": (
                    None if level1_manifest is None else level1_manifest.get("mu_sha256")
                ),
                "n_estimated": None if level1_params is None else level1_params.n_estimated,
                "n_unestimated": None if level1_params is None else level1_params.n_unestimated,
                "sigma_min": level1_sigma_min if include_robust_z else None,
            }
            if include_robust_z
            else {"enabled": False}
        ),
    }
    # ADR 0008: orient gene-panel MBS (exclude residual slot) vs signed gene-mean M.
    polarity = "hyper_aligned"
    if train_records:
        model.eval()
        head.eval()
        mbs_rows: list[np.ndarray] = []
        present_rows: list[np.ndarray] = []
        m_batches: list[np.ndarray] = []
        gene_batches: list[np.ndarray] = []
        for rec in train_records:
            batch = pack_hier_records_to_batch(
                [rec],
                locus_region=locus_region,
                age_values=[None],
                age_enabled=[False],
                tissue_enabled=[True],
            ).to(device)
            with torch.no_grad():
                mbs, present = _packed_hier_mbs(model, batch)
            gene_mbs = mbs[:, :n_genes].detach().cpu().numpy()[0]
            gene_present = present[:, :n_genes].detach().cpu().numpy()[0]
            mbs_rows.append(gene_mbs)
            present_rows.append(gene_present)
            feats = rec.features.cpg_features
            regions = rec.features.cpg_to_region
            # flat_standard layout: beta col 0, M col 1 (see feature_schema.FLAT_STANDARD).
            if feats.shape[0] and feats.shape[1] > 1 and regions.size:
                m_batches.append(np.asarray(feats[:, 1], dtype=np.float64))
                gene_batches.append(
                    np.asarray(
                        locus_region.region_to_gene[regions.astype(np.int64)],
                        dtype=np.int64,
                    )
                )
        if mbs_rows and m_batches:
            oriented = orient_run_scores(
                np.stack(mbs_rows, axis=0),
                signed_m=accumulate_signed_gene_mean_m(
                    n_genes=n_genes,
                    cpg_m_batches=m_batches,
                    cpg_to_gene_batches=gene_batches,
                ),
                present=np.stack(present_rows, axis=0),
            )
            polarity = str(oriented["score_polarity"])
    manifest = score_manifest(score_polarity=polarity, fold_id=None, restart_id=run_id)
    (run_root / "score_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    metrics_out["score_manifest"] = manifest
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
        "n_panel": n_panel,
        "n_regions": locus_region.n_regions,
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
        json.dumps(metrics_out, indent=2, default=str) + "\n",
        encoding="utf-8",
    )

    tb_url = None
    tb_port_out = None
    hint = None
    if tb_server is not None:
        tb_url = tb_server.url
        tb_port_out = tb_server.port
        hint = f"TensorBoard: {tb_url}"

    return TrainResult(
        run_id=run_id,
        run_dir=run_root,
        checkpoint_dir=ckpt_root,
        metrics=metrics_out,
        best_epoch=best_epoch,
        tensorboard_url=tb_url,
        tensorboard_port=tb_port_out,
        monitor_hint=hint,
    )
