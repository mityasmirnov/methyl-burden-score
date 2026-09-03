# Plan: 7G′ Stage A DeepRVAT-style architecture screen

Status: **in progress** — Stage A reopened for compute-efficient screening
(mixed pooling, RBS diagnostic, vector cascade, one-hop). Required P2/P4/P5-max
/`C-mvalue-*-G` GPU runs already landed; cascade is **not** ≥0.03 ahead of
classical. **P5 inactive.** Stage B GPU and Milestone **7** remain blocked.

Last updated: **2026-09-03** (orientation contract v2 fix; L1/L5 training configs;
`early_stopping_start_epoch`; annotation ablation grid).

Parent: [`milestone-7g-prime-matched-probe-lightweight.md`](milestone-7g-prime-matched-probe-lightweight.md).
Normative: [ADR 0010](../adr/0010-gene-allocation-policy.md),
[`ARCHITECTURE.md`](../ARCHITECTURE.md), [`SCORING_PIPELINE.md`](../SCORING_PIPELINE.md).

## Scope and acceptance

Expand Stage A into a small architecture-screening benchmark on the existing
`explicit_only` gene-linked panel (51 375 CpGs) and frozen
`hub-ats-7e-3fold-v1` folds. Always report tissue, age, and sex. Select by
per-task rankings + Pareto, not tissue F1 alone.

**Done when:** promoted arms have three valid held-out folds; the report
identifies where age/sex information is lost; one gene-aggregation architecture
is selected **or** the report concludes no neural aggregation is yet adequate.

## Locked decisions

| Choice | Decision | Why |
|--------|----------|-----|
| Panel | Reuse `explicit_only` + `gene_panel_manifest.json` | Fair vs classical `-G` |
| P5 | Inactive; keep historical artifacts | 30 epochs hurt; no more epoch grid |
| Scalar pooling | Independent `cpg_pool` / `region_pool` | P4 changed both levels at once |
| Vector cascade | Pool region **embeddings**, then `rho_G` | Scalar RBS before gene pool discards promoter/body |
| One-hop | Annotated CpG → gene pool (batched flat train) | Closest DeepRVAT literal copy |
| RBS diagnostic | `all_gene_rbs.zarr` (not orphan `rbs.zarr`) | Locate age/sex loss before vs after gene pool |
| Screening | Tier 1: 5 ep; Tier 2: 15 ep if Pareto/near-best | Compute-efficient |
| Execution | **Model-by-model**; refresh `analysis.md` after each arm | Start with one-hop (`N-light-gene-*`), then cascade screen |
| Regulatory cCRE | Slots reserved; zeros until graph release | SCREEN deferred (7C non-goal) |
| Stage B seed-gene | Design only; separate from fold-panel Stage B | Do not delay Stage A |
| Orientation contract | **v2** (2026-09-03): `evaluate_flat_mbs_e2e` passes **raw** encoder MBS to heads; `orient_mbs_array` affects only the exported association artifact. Legacy checkpoints with negated head weights use the repair path (`legacy_negated_heads=True`). | Passing `1-MBS` through unchanged head weights gives wrong logits (ADR 0008). |
| Training LR baseline | **L1**: single `1e-3` LR for all parameters (no `head_lr_multiplier`) — DeepRVAT-like. **L5**: head LR 5× diagnostic, only if L1 representation diagnostics justify it. | Avoids premature optimisation decisions before representation quality is understood. |
| Early stopping | `early_stopping_start_epoch: 5`; patience 3–5 | Short warmup prevents premature convergence signal on first epochs. |
| Annotation channels | `cpg_context` wired (UCSC CGI from `loci.parquet`). Gene-role one-hot populated. Regulatory (cCRE/DHS/ChromHMM) reserved zero slots — sources not on disk. | Stage A non-goal for regulatory; ablation grid tests what is available. |

## Arm naming

| Arm | Role |
|-----|------|
| `N-cascade-scalar-max-max` | Alias of landed `P2-G` |
| `N-cascade-scalar-mean-mean` | Alias of landed `P4-G` |
| `N-cascade-scalar-mean-max` | New mixed train |
| `N-cascade-scalar-max-mean` | New mixed train |
| `N-cascade-vector-mean-max` | Preferred vector hypothesis |
| `N-cascade-vector-max-max` | Vector max/max |
| `N-light-gene-max` / `N-light-gene-mean` | One-hop annotated DeepSet |
| `rbs_linear_probe` / `rbs_enet` | Frozen gene-linked RBS readouts |

## Data / artifact flow

```mermaid
flowchart LR
  panel["gene_panel_manifest.json"] --> scalar["scalar cascade 2x2"]
  panel --> vector["vector cascade"]
  panel --> light["one-hop FlatDeepSetRegion"]
  panel --> classical["C-mvalue-*-G"]
  scalar --> rbsDiag["all_gene_rbs.zarr"]
  scalar --> report["stage0_7g_gene_only_probe"]
  vector --> report
  light --> report
  classical --> report
  rbsDiag --> report
```

## Non-goals

- Milestone **7** 5×6 OOF.
- Additional P5 / 30-epoch arms.
- SCREEN/cCRE graph rebuild.
- Joint end-to-end orphan/direct fusion in Stage A.
- Implementing Stage B seed-gene transfer training (design only).

## Stage B seed-gene transfer (design only — do not implement training now)

DeepRVAT discovers trait-associated **seed genes**, trains a **shared** impairment
function through phenotype heads on those genes, cross-fits **across samples**,
then applies the shared function to seed **and** non-seed genes, and evaluates
association discovery/replication. It is **not** ordinary train-gene / test-gene CV.

Keep the existing fold-selected-panel Stage B (`stage0_7g_prime_stage_b.yaml`,
`C-mvalue-enetS` / `N-cascade-S`) as a **separate** question. Probe selection and
gene-transfer validation must not be conflated.

### Protocol sketch

1. **Inspect coverage** — Hub phenotype labels + EWAS Atlas CpG→trait tables;
   record sample counts, study IDs, and gene mappings under
   `reports/inspection/` (do not invent a trait).
2. **Trait bar** — propose **one** trait with: adequate labeled samples on the
   Hub matrix; a substantial Atlas CpG set; enough distinct genes after
   `explicit_only` mapping; ≥2 independent studies for replication. If none
   meet the bar, stop and document.
3. **Seeds** — significant Atlas CpGs with defensible gene edges (same ADR 0010
   policy as Stage A). Persist seed gene IDs + provenance.
4. **Study overlap control** — exclude or flag seed evidence from studies in the
   evaluation fold; prefer independent seed vs validation study sets.
5. **Train** — shared encoder on **seed genes only**; sample-level cross-fitting
   unchanged (gene transfer does **not** replace held-out samples).
6. **Apply** — score all eligible genes with the shared function.
7. **Evaluate non-seed transfer** — enrichment for held-out Atlas associations;
   replication in independent studies; permutation calibration; vs regional
   means and classical CpG burdens.
8. **Robustness** — drop 10% of seeds; add 20% non-associated genes; vary seed
   count.

### Atlas / Hub coverage note (inspect only; no trait invented)

From `reports/inspection/deepmat_data_v1/trait_eligibility.md` (2026-09-03):

- Hub ATS core phenotypes (age / tissue / sex) are already Stage A targets — **not**
  seed-gene transfer traits.
- Pack-level `disease` / `cancer` binary labels still fail the project bar
  (need ≥200 cases, ≥200 controls, ≥3 studies with unknown≠control).
- Candidate **continuous** packs with multi-study support include BMI
  (`bmi` pack: n≈2070, 25 studies) and blood/brain tissue packs — but seed-gene
  transfer needs an **EWAS Atlas CpG set** for that trait, not just Hub labels.
- **Next inspect before proposing a trait:** join Atlas association tables to Hub
  sample study IDs; count distinct genes after `explicit_only`; check independent
  discovery vs replication studies. If no trait clears that join, document stop.

### Non-goals for this design note

- Starting Milestone **7** OOF.
- Replacing Stage B fold-selected panel GPU work.
- Fetching SCREEN/cCRE for regulatory features (still reserved zeros).

## N-light one-hop status and repair (2026-09-03)

### What was wrong

Pre-fix `N-light-gene-*` `mbs_e2e` numbers (~0.000–0.001 tissue F1) were
invalid due to two compounding bugs:

1. **Orientation anchor**: `orient_mbs_array` was called **before** the
   phenotype heads, passing `1 − MBS` through head weights that expected raw
   scores → wrong logits at evaluation.
2. **Head/score mismatch**: `_orient_and_write_score_manifest()` was also
   mutating head weights at checkpoint time, creating inconsistent state.

Both are fixed in commit `fc8cd6f` (2026-09-03). Frozen MBS probes (linear /
elastic-net on raw MBS) showed valid representation signal throughout — one-hop
is **not** rejected on the basis of the broken `mbs_e2e` numbers.

### Corrected evaluation contract (v2)

```
evaluate_flat_mbs_e2e:
  heads receive: raw encoder MBS (always)
  exported association artifact: orient_mbs_array(mbs, score_polarity)
  legacy repair: legacy_negated_heads=True → 1-MBS (logit-preserving for that ckpt)
```

Unit tests: `tests/unit/test_orientation_eval.py` —
`test_v2_contract_heads_see_raw_mbs_when_flipped`,
`test_legacy_path_heads_see_one_minus_mbs_when_flipped`,
`test_hyper_aligned_heads_always_see_raw`.

### GPU run sequence for one-hop (fold 0)

Run in this order; do not start annotation ablations before L1 numbers are in.

1. **L1 baseline** (`light_mean_l1.yaml`) — M-only features, mean pooling,
   single LR `1e-3`, `early_stopping_start_epoch: 5`, patience 3, max 10 ep.
   Establishes the DeepRVAT-like single-LR control.
2. **Inspect representation diagnostics** — gene MBS variance, head weight
   norms, validation curves. If encoder variance is reasonable, proceed to
   ablations. If heads are flat despite encoder signal, run L5.
3. **L5 head-LR diagnostic** (`light_mean_l5.yaml`) — identical to L1 but
   `head_lr_multiplier: 5.0`. Run only if L1 diagnostics justify it.
4. **Annotation ablation grid** (all 18 arms, `run_7g_gene_only_probe.py
   --fold 0 --device cuda`):

| Arm | Config suffix | Feature mode | Notes |
|-----|--------------|-------------|-------|
| A0 | `ablation_m_only` | `m_only` | M + observed; no annotation |
| A1 | `ablation_m_role` | `m_role` | M + gene-role one-hot |
| A2 | `ablation_m_context` | `m_context` | M + CpG-context one-hot |
| A3 | `ablation_m_role_context` | `m_role_context` | M + role + context |
| A4/A7 | `ablation_full` | `full` | All channels (reg stays zero) |
| N0 | `ablation_n0_obs_only` | `obs_only` | Observed flag only; no M |
| N1 | `ablation_n1_anno_only` | `anno_only` | Role + context; no M |
| N2 | `ablation_n2_reg_permuted` | `full` + permuted reg | Regulatory permutation sanity |
| N3 | `ablation_n3_reg_zero` | `full` + zero reg | Explicit zero-reg control |

Each arm has two seeds (`_s2` suffix). Report bootstrap CIs and representation
diagnostics via `write_7g_gene_only_probe_report.py`.

5. **Post-hoc direct/orphan fusion** — deferred until ablation grid is scored.

## Open questions (resolved by this screen)

1. Mean vs max at each cascade level.
2. Whether scalar RBS discards information.
3. Whether gene pooling discards information relative to RBS.
4. Whether one-hop matches or beats the cascade.
5. Whether one scalar MBS per gene is enough for age/sex/tissue.
6. Whether `cpg_context` annotation adds measurable signal over M-only baseline (A0 vs A2/A3). ← *new, added 2026-09-03*
