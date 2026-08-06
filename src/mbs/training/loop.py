"""Plain PyTorch training loop for the flat DeepRVAT-style baseline."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import torch
import torch.nn.functional as F
import yaml
from torch import nn

from mbs.evaluation.splits import build_study_grouped_split
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
    record_to_batch,
)
from mbs.training.features import build_static_column_table
from mbs.training.locus_gene import LocusGeneIndex, build_locus_gene_index, load_graph_tables
from mbs.training.phenotypes import SamplePhenotype, load_gse35069_phenotypes
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


def _forward_tissue_loss(
    *,
    model: FlatDeepSet,
    head: SeedMaskedLinearHead,
    batch: FlatBatch,
    class_weights: torch.Tensor | None,
) -> tuple[torch.Tensor, dict[str, float]]:
    output = model(batch.cpg_features, batch.cpg_to_gene, batch.n_genes)
    mbs = output["mbs"].unsqueeze(0)
    present = output["present"].unsqueeze(0)
    logits = head(mbs, present)
    tissue_loss = F.cross_entropy(
        logits,
        batch.tissue_target,
        weight=class_weights.to(logits.device) if class_weights is not None else None,
    )
    pred = int(logits.argmax(dim=-1).item())
    metrics = {
        "loss": float(tissue_loss.detach().item()),
        "tissue_correct": float(pred == int(batch.tissue_target.item())),
    }
    return tissue_loss, metrics


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


def _run_epoch(
    *,
    records: list[FlatSampleRecord] | None,
    phenotypes: list[SamplePhenotype] | None,
    pilot_store: _PilotStore | None,
    model: FlatDeepSet,
    head: SeedMaskedLinearHead,
    optimizer: torch.optim.Optimizer | None,
    device: torch.device,
    n_genes: int,
    class_weights: torch.Tensor | None,
    use_amp: bool,
    grad_clip: float,
    train: bool,
) -> dict[str, float]:
    if train:
        model.train()
        head.train()
    else:
        model.eval()
        head.eval()

    total_loss = 0.0
    total_correct = 0.0
    n = 0
    amp_dtype = torch.bfloat16 if device.type == "cuda" else torch.float32
    context = torch.enable_grad() if train else torch.no_grad()

    def _iter_records() -> list[FlatSampleRecord]:
        if records is not None:
            return records
        if phenotypes is None or pilot_store is None:
            raise RuntimeError("pilot phenotypes and store are required when records is None")
        return [_materialize_record(ph, pilot_store) for ph in phenotypes]

    with context:
        for record in _iter_records():
            batch = record_to_batch(record, n_genes=n_genes).to(device)
            if train and optimizer is not None:
                optimizer.zero_grad(set_to_none=True)
            with torch.autocast(device_type=device.type, dtype=amp_dtype, enabled=use_amp):
                loss, metrics = _forward_tissue_loss(
                    model=model,
                    head=head,
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
            total_loss += metrics["loss"]
            total_correct += metrics["tissue_correct"]
            n += 1

    return {
        "loss": total_loss / max(n, 1),
        "accuracy": total_correct / max(n, 1),
        "n_samples": float(n),
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
    model_cfg = config.get("model", {})

    fixture_records: list[FlatSampleRecord] | None = None
    train_phenotypes: list[SamplePhenotype] | None = None
    val_phenotypes: list[SamplePhenotype] | None = None
    pilot_store: _PilotStore | None = None
    train_records: list[FlatSampleRecord] | None = None
    val_records: list[FlatSampleRecord] | None = None
    class_weights: torch.Tensor | None = None

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
    else:
        pilot = config.get("pilot", {})
        matrix_id = str(pilot["matrix_id"])
        graph_id = str(pilot["graph_id"])
        feature_set = str(
            pilot.get("static_feature_set") or config.get("features", {}).get("static_feature_set")
        )
        meta_rel = Path(str(pilot["phenotype_metadata"]))
        metadata_path = meta_rel if meta_rel.is_absolute() else data_root / meta_rel

        matrix_paths = matrix_store_paths(data_root / "canonical" / "matrices" / matrix_id)
        sample_index = read_sample_index(matrix_paths.sample_index_path)
        locus_index = read_locus_index(matrix_paths.locus_index_path)
        sample_ids = sample_index.sort_values("row_index")["sample_id"].astype(str).tolist()
        phenotypes, class_names = load_gse35069_phenotypes(metadata_path, sample_ids=sample_ids)

        train_donors = {str(x) for x in pilot.get("train_donors", ["1", "2", "3", "4"])}
        val_donors = {str(x) for x in pilot.get("val_donors", ["5", "6"])}
        train_phenotypes, val_phenotypes, split = _split_by_donor(
            phenotypes, train_donors=train_donors, val_donors=val_donors
        )

        lr_edges, regions = load_graph_tables(data_root / "canonical" / "graphs" / graph_id)
        locus_gene = build_locus_gene_index(
            locus_index=locus_index,
            locus_region_edges=lr_edges,
            regions=regions,
            max_loci=max_loci,
        )
        gene_ids = locus_gene.gene_ids
        n_genes = locus_gene.n_genes
        n_classes = len(class_names)
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
        class_weights = _class_weights([p.class_index for p in train_phenotypes], n_classes)
        if torch.allclose(class_weights, torch.ones_like(class_weights)):
            class_weights = None
        split.update(
            {
                "n_genes": n_genes,
                "n_classes": n_classes,
                "class_names": class_names,
                "matrix_id": matrix_id,
                "max_loci": max_loci,
            }
        )

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
    seed_mask = torch.ones(n_classes, n_genes, dtype=torch.float32, device=device)
    head = SeedMaskedLinearHead(n_genes, n_classes, seed_mask).to(device)
    # Zero gene_weight is a saddle for CE on saturated MBS; break symmetry.
    with torch.no_grad():
        head.gene_weight.normal_(0.0, 0.05)
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
    tb_writer = None
    jsonl_path = run_root / "metrics.jsonl"
    if use_tb:
        from torch.utils.tensorboard import SummaryWriter  # noqa: PLC0415

        tb_dir = run_root / "tb"
        tb_dir.mkdir(parents=True, exist_ok=True)
        tb_writer = SummaryWriter(log_dir=str(tb_dir))

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
        )
        row = {
            "epoch": epoch,
            "train_loss": train_metrics["loss"],
            "train_accuracy": train_metrics["accuracy"],
            "val_loss": val_metrics["loss"],
            "val_accuracy": val_metrics["accuracy"],
            "learning_rate": lr,
        }
        history.append(row)
        with jsonl_path.open("a", encoding="utf-8") as jsonl_handle:
            jsonl_handle.write(json.dumps(row) + "\n")
        if tb_writer is not None:
            tb_writer.add_scalar("loss/train", row["train_loss"], epoch)
            tb_writer.add_scalar("loss/val", row["val_loss"], epoch)
            tb_writer.add_scalar("accuracy/train", row["train_accuracy"], epoch)
            tb_writer.add_scalar("accuracy/val", row["val_accuracy"], epoch)
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
    )
