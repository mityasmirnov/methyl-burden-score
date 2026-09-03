"""Canonical labels for Stage 0 evaluation arm IDs (7E–7G′ reports)."""

from __future__ import annotations

from typing import Iterable

# Normative reference: docs/ARM_GLOSSARY.md (generated from this module).

ARM_DESCRIPTIONS: dict[str, str] = {
    # --- 7G tissue probe (P = probe ablation on cascade; test metric often late fusion) ---
    "P0-baseline": (
        "Replay 7G product cascade (`N-cascade-l1`): train + **late fusion** "
        "`[orphan_rbs | mbs | direct_contrib]` on all 65 536 prefix loci."
    ),
    "P1-fusion-tissue-heavy": (
        "**Fusion-only** on saved P0 scores: balanced logistic tissue head + PCA(32); "
        "no retrain."
    ),
    "P2-end2end-tissue-weight": (
        "Cascade **retrain** with tissue_loss_weight=3, age_loss_weight=0.3; "
        "test via **late fusion** (not MBS-only)."
    ),
    "P2-fusion-balanced": "P2 saved scores + balanced logistic late fusion (no retrain).",
    "P2-fusion-logistic": "P2 saved scores + multinomial logistic late fusion (no retrain).",
    "P2-fusion-sgd": "P2 saved scores + SGD logistic late fusion (no retrain).",
    "P3-region-head-bypass": (
        "Transparent **`T-mean-region`**: region-mean M-values → linear heads; "
        "no CascadeDeepSet training."
    ),
    "P4-pooling-mean": (
        "Like P2 loss weights but **mean/mean** CpG→region and region→gene pooling; "
        "late-fusion test."
    ),
    "P4-fusion-balanced": "P4 saved scores + balanced logistic late fusion.",
    "P4-fusion-logistic": "P4 saved scores + logistic late fusion.",
    "P4-fusion-sgd": "P4 saved scores + SGD late fusion.",
    "P5-epochs-30": (
        "P2 loss weights, **max/max** pooling, **30-epoch** ceiling + early stop on "
        "validation tissue F1; late-fusion test."
    ),
    "P5-mean-epochs-30": "P5 schedule with **mean/mean** pooling.",
    "P5-fusion-balanced": "P5 saved scores + balanced logistic late fusion.",
    # --- 7G′ Stage A (suffix -G = gene-linked CpG panel only; primary metric mbs_e2e) ---
    "P2-G": (
        "CascadeDeepSet on **gene-linked CpGs only** (`gene_cols`); max/max pooling; "
        "15 epochs; P2 loss weights; primary metric **`mbs_e2e`** (end-to-end MBS heads)."
    ),
    "P4-G": (
        "Like P2-G but **mean/mean** pooling; same gene-only panel and **`mbs_e2e`** metric."
    ),
    "P5-G-max": (
        "Gene-only cascade; max/max pooling; 30-epoch cap + early stop; **`mbs_e2e`** metric. "
        "**Inactive** for Stage A DeepRVAT screen (historical only)."
    ),
    "P5-G-mean": (
        "Gene-only cascade; mean/mean pooling; 30-epoch cap + early stop; **`mbs_e2e`** metric. "
        "**Inactive**."
    ),
    "N-cascade-scalar-max-max": "Alias of **`P2-G`**: scalar RBS cascade, CpG max → gene max.",
    "N-cascade-scalar-mean-mean": "Alias of **`P4-G`**: scalar RBS cascade, CpG mean → gene mean.",
    "N-cascade-scalar-mean-max": (
        "Scalar RBS cascade; **mean** CpG→region + **max** region→gene; Stage A screen."
    ),
    "N-cascade-scalar-max-mean": (
        "Scalar RBS cascade; **max** CpG→region + **mean** region→gene; Stage A screen."
    ),
    "N-cascade-vector-mean-max": (
        "Vector-region cascade: pool region **embeddings** (mean CpG→region, max region→gene) "
        "then ``rho_G`` → MBS."
    ),
    "N-cascade-vector-max-max": (
        "Vector-region cascade with max/max pooling of latent region embeddings."
    ),
    "N-light-gene-max": (
        "One-hop annotated FlatDeepSetRegion: `[M, gene-role, CGI context, regulatory, flags]` "
        "→ gene **max** → MBS."
    ),
    "N-light-gene-mean": (
        "One-hop annotated FlatDeepSetRegion with gene **mean** pooling."
    ),
    "P2-orphan-ablation": (
        "Full assignment (not gene-only): same P2 loss weights; compares "
        "**`fusion_full`** (orphan RBS + MBS + direct) vs **`fusion_mbs_direct`** "
        "(MBS + direct, orphan block dropped)."
    ),
    "C-mvalue-classical-G": (
        "Bundle runner for all **`C-mvalue-*-G`** arms below on the same `gene_cols` panel."
    ),
    # --- Classical M-value baselines (C = classical; methylation input only) ---
    "C-mvalue-ridge": (
        "All prefix loci → M-value → train-fold impute/scale → ridge (age) / "
        "SGD-L2 logistic (tissue, sex)."
    ),
    "C-mvalue-ridge-G": "Same as `C-mvalue-ridge` but only **`gene_cols`** (56k matched panel).",
    "C-mvalue-enet": (
        "All prefix loci → M-value → elastic-net regression (age) / logistic enet "
        "(tissue, sex); **7G bake-off tissue winner** (F1≈0.334 on 65k)."
    ),
    "C-mvalue-enet-G": (
        "Elastic-net on **identical gene-linked CpGs** as neural `-G` arms (fair Stage A comparator)."
    ),
    "C-mvalue-enetS": (
        "**Stage B:** fold-safe stability-selected sparse CpG panel (outer-train only) → "
        "refit enet; one test eval per fold."
    ),
    "C-mvalue-hgb": "M-value → HistGradientBoosting (stands in for LightGBM) per trait.",
    "C-mvalue-hgb-G": "HGB on **`gene_cols`** only.",
    "C-mvalue-sva": "M-value → PCA surrogate components (Bioconductor `sva` stand-in) → ridge/logistic.",
    "C-mvalue-sva-G": "PCA-SVA + ridge on **`gene_cols`** only.",
    # --- Neural product / Stage B (N = neural encoder) ---
    "N-cascade-l1": (
        "**7F/7G product topology:** CascadeDeepSet (CpG→RBS→MBS) + fold-fitted direct enet "
        "+ **late fusion** on saved blocks; Level-1 z on direct branch."
    ),
    "N-cascade-S": (
        "**Stage B:** locked Stage-A cascade hyperparameters on the **fold-selected** CpG panel "
        "(same loci as `C-mvalue-enetS` per fold)."
    ),
    "N-light-type": (
        "Legacy Stage B name for FlatDeepSetRegion; prefer **`N-light-gene-max`** / "
        "**`N-light-gene-mean`** in Stage A screen."
    ),
    "N-mbs-posthoc-full-fusion": (
        "**Stage B post-hoc fusion:** MBS encoder trained once; CPU late fusion on "
        "**`[orphan_rbs | mbs | direct_contrib]`** (`fusion_full`). Not joint end-to-end."
    ),
    "N-mbs-posthoc-mbs-direct": (
        "**Stage B orphan ablation:** same encoder; CPU fusion on **`[mbs | direct_contrib]`** only "
        "(`fusion_mbs_direct`)."
    ),
    "N-full": (
        "Deprecated alias for **`N-mbs-posthoc-full-fusion`** (post-hoc fusion, not joint training)."
    ),
    "N-mbs-direct-only": (
        "Deprecated alias for **`N-mbs-posthoc-mbs-direct`**."
    ),
    # --- Transparent linear baselines (T = transparent; fixed features, no shared encoder train) ---
    "T-mean-gene": "Presence-aware **gene-mean** M-values → fold-fitted linear multitask heads.",
    "T-mean-region": "Presence-aware **region-mean** M-values → linear heads (strong tissue baseline).",
    "T-enet": "Transparent elastic-net on gene- or region-aggregated M-value features.",
    # --- 7E development CV arm ids (historical) ---
    "N-multipath-l1a": (
        "7E budget-limited multipath arm (region-mean linear fusion; **not** 7F neural product)."
    ),
    "N-flat-l1": "7E flat DeepSet + Level-1 normalization variant.",
    "N-hier-l1": "7E hierarchical DeepSet + Level-1 variant.",
}

EVAL_MODE_DESCRIPTIONS: dict[str, str] = {
    "mbs_e2e": (
        "End-to-end **`MultitaskHeads`** on MBS only (no orphan/direct columns) — "
        "**Stage A primary metric**."
    ),
    "mbs_linear_probe": (
        "CPU linear probe fit on **saved MBS matrix** only (representation check, not fusion)."
    ),
    "mbs_enet": (
        "CPU elastic-net heads on the **same saved MBS** (age + tissue + sex). "
        "Same CascadeDeepSet weights as `mbs_e2e`; no encoder retrain."
    ),
    "rbs_linear_probe": (
        "CPU linear probe on frozen **gene-linked RBS** (`all_gene_rbs.zarr`); "
        "diagnoses loss before vs after gene pooling."
    ),
    "rbs_enet": (
        "CPU elastic-net on frozen gene-linked RBS (same matrix as `rbs_linear_probe`)."
    ),
    "fusion_full": "Late fusion on **`[orphan_rbs | mbs | direct_contrib]`** columns.",
    "fusion_mbs_direct": "Late fusion on **`[mbs | direct_contrib]`** — orphan RBS ablation.",
}


def arm_description(arm_id: str) -> str:
    """Return human-readable description for an arm id (or a short fallback)."""
    if arm_id in ARM_DESCRIPTIONS:
        return ARM_DESCRIPTIONS[arm_id]
    return f"See docs/ARM_GLOSSARY.md (no inline entry for `{arm_id}`)."


def render_prefix_legend() -> list[str]:
    """Explain arm-id prefixes used across 7E–7G′ reports."""
    return [
        "### Naming prefixes",
        "",
        "| Prefix | Meaning | Examples |",
        "|--------|---------|----------|",
        "| **P** | **Probe** ablation (7G tissue investigation or 7G′ `-G` gene-only) | `P2-end2end-tissue-weight`, `P4-G` |",
        "| **C-mvalue-** | **Classical** sklearn/HGB on M-values (no neural encoder) | `C-mvalue-enet`, `C-mvalue-enet-G` |",
        "| **N-** | **Neural** encoder (CascadeDeepSet or flat/hier variants) | `N-cascade-l1`, `N-light-type` |",
        "| **T-** | **Transparent** fixed-feature linear baselines | `T-mean-region`, `T-mean-gene` |",
        "| **-G** | **Gene-linked** CpG panel only (`gene_cols`; 7G′ Stage A) | `P2-G`, `C-mvalue-ridge-G` |",
        "| **-S** | Fold-safe **Selected** sparse panel (7G′ Stage B) | `C-mvalue-enetS`, `N-cascade-S` |",
        "",
    ]


def render_arm_glossary_section(
    arm_ids: Iterable[str],
    *,
    title: str = "## Arm glossary",
    include_prefix_legend: bool = True,
    extra_eval_modes: Iterable[str] | None = None,
) -> list[str]:
    """Markdown lines: legend + table of arm id → description for arms in this report."""
    ids = list(dict.fromkeys(arm_ids))  # preserve order, unique
    lines: list[str] = [title, ""]
    if include_prefix_legend:
        lines.extend(render_prefix_legend())
    lines.extend(
        [
            "| Arm | Description |",
            "|-----|-------------|",
        ]
    )
    for arm_id in ids:
        desc = ARM_DESCRIPTIONS.get(arm_id, arm_description(arm_id))
        lines.append(f"| `{arm_id}` | {desc} |")
    if extra_eval_modes:
        modes = list(dict.fromkeys(extra_eval_modes))
        if modes:
            lines.extend(["", "### Evaluation modes (cascade metrics)", ""])
            for mode in modes:
                lines.append(f"- **`{mode}`**: {EVAL_MODE_DESCRIPTIONS.get(mode, mode)}")
    lines.append("")
    return lines


def render_full_glossary_doc() -> str:
    """Standalone markdown reference for all known arm ids and eval modes."""
    lines = [
        "# Evaluation arm glossary (7E–7G′)",
        "",
        "Canonical source: `src/mbs/inspection/arm_glossary.py`. "
        "Inspection reports embed a subset via `render_arm_glossary_section()`.",
        "",
    ]
    lines.extend(render_prefix_legend())
    lines.extend(["## All arms", "", "| Arm | Description |", "|-----|-------------|"])
    for arm_id in sorted(ARM_DESCRIPTIONS):
        lines.append(f"| `{arm_id}` | {ARM_DESCRIPTIONS[arm_id]} |")
    lines.extend(["", "## Evaluation modes (cascade metrics)", ""])
    for mode, desc in sorted(EVAL_MODE_DESCRIPTIONS.items()):
        lines.append(f"- **`{mode}`**: {desc}")
    lines.append("")
    return "\n".join(lines)
