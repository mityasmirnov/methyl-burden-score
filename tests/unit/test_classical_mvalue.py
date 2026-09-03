"""Classical M-value fold eval: elastic-net age must stay in year scale."""

from __future__ import annotations

import numpy as np

from mbs.training.classical_mvalue import fit_eval_mvalue_fold


def test_enet_age_mae_not_blanked_on_year_scale_high_dim() -> None:
    rng = np.random.default_rng(0)
    n_tr, n_te, p = 200, 80, 4000
    weights = np.zeros(p)
    weights[:8] = rng.normal(size=8)
    x_tr = rng.normal(size=(n_tr, p)).astype(np.float32)
    x_te = rng.normal(size=(n_te, p)).astype(np.float32)
    age_tr = 40.0 + x_tr @ weights + rng.normal(scale=3.0, size=n_tr)
    age_te = 40.0 + x_te @ weights + rng.normal(scale=3.0, size=n_te)
    empty_int = np.zeros(n_tr, dtype=np.int64)
    ph_tr = {
        "age": age_tr,
        "age_mask": np.ones(n_tr, dtype=bool),
        "tissue": empty_int,
        "tissue_mask": np.zeros(n_tr, dtype=bool),
        "sex": empty_int,
        "sex_mask": np.zeros(n_tr, dtype=bool),
    }
    ph_te = {
        "age": age_te,
        "age_mask": np.ones(n_te, dtype=bool),
        "tissue": np.zeros(n_te, dtype=np.int64),
        "tissue_mask": np.zeros(n_te, dtype=bool),
        "sex": np.zeros(n_te, dtype=np.int64),
        "sex_mask": np.zeros(n_te, dtype=bool),
        "tissues": np.array(["x"] * n_te),
    }
    out = fit_eval_mvalue_fold(x_tr, x_te, ph_tr, ph_te, "enet", impute=False)
    assert out.get("age") is not None, out.get("age_note")
    mae = float(out["age"]["mae"])
    assert mae < 100.0
    assert np.isfinite(mae)
