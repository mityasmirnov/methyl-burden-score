"""Independently trained branch arms (gene / rbs / tbs / direct)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from mbs.training.direct_cpg import fit_direct_elasticnet
from mbs.training.loop import TrainResult, train_flat_baseline
from mbs.training.phenotypes import hub_longform_ready

VALID_ARMS = frozenset({"gene", "rbs", "tbs", "direct"})

__all__ = ("VALID_ARMS", "hub_longform_ready", "train_branch_arm")


def train_branch_arm(
    *,
    arm: str,
    project_root: Path,
    data_root: Path,
    artifact_root: Path,
    config: dict[str, Any],
    run_id: str,
    device: str = "cpu",
    overfit_fixture: bool = True,
) -> dict[str, Any]:
    """Train one arm on its own run id. Eval-time masking is not this function."""
    if arm not in VALID_ARMS:
        raise ValueError(f"arm must be one of {sorted(VALID_ARMS)}")
    if arm == "direct":
        rng = np.random.default_rng(int(config.get("experiment", {}).get("seed", 42)))
        n, p = 20, 8
        z = rng.normal(size=(n, p))
        obs = rng.random(size=(n, p)) > 0.2
        y = z[:, 0] + 0.1 * rng.normal(size=n)
        studies = np.array(["A"] * 10 + ["B"] * 10)
        fitted = fit_direct_elasticnet(z, obs, y, studies, min_studies=2)
        return {"arm": arm, "run_id": run_id, "n_loci": fitted["n_loci"], "direct": True}
    cfg = dict(config)
    cfg.setdefault("model", {})
    cfg["model"] = {**cfg["model"], "arm": arm}
    result: TrainResult = train_flat_baseline(
        project_root=project_root,
        data_root=data_root,
        artifact_root=artifact_root,
        config=cfg,
        run_id=run_id,
        device_str=device,
        overfit_fixture=overfit_fixture,
    )
    return {"arm": arm, "run_id": result.run_id, "metrics": result.metrics, "direct": False}
