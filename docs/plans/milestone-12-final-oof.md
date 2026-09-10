# Milestone 12: Final study-grouped OOF (5×6)

> **Alias:** historical **Milestone 7** in ADR 0007 / older docs.
> Index: [`MILESTONE_INDEX.md`](MILESTONE_INDEX.md).
>
> **Renumber note (2026-09-07):** previously stubbed as Milestone **13**;
> expression auxiliary is now **13** (after OOF).

**Status:** **N-light arm `in_progress` ← NOW** — GPU 2 at **30 epochs** +
n>200 aux; **24/30** complete (folds 0–3); **f4** in flight. 16-ep plumbing
archived (`stage0-12-nlight-oof-ep16-f0-r0`; nested **0.300 / 8.63 / 0.880**).

**Question:** Under nested enet, does N-light@64 hold on study-grouped 5×6 at
the **65k-prefix** HM450 panel?

**Approaches tested:** `stage0_12_nlight_oof.yaml` / `run_12_nlight_oof.sh`;
30-ep + n>200 aux; product readout `mbs_enet_nested`.

**Results (interim 24/30):** nested **0.325 / 9.79 / 0.831** (tissue / age MAE
/ sex AUROC); e2e **0.319 / 15.17 / 0.814**. Vs 3-fold nested ref
**0.368 / 9.88 / 0.803**. Report:
[`../../reports/inspection/stage0_12_nlight_oof/analysis.md`](../../reports/inspection/stage0_12_nlight_oof/analysis.md).
Final numbers `pending` until `check_12_oof_completeness.py` is 30/30.

**Verdict:** open — interim on-track; 65k **validation** only.

**Cascade arm:** **blocked on post–N-light GATE G1–G4** (gene utilization,
positional/CpGPT, platform robustness, fair 10e + 10d). Default topology
**native P2-G** unless fair **10e** flips it. No auto 65k-matched cascade
queue. Does **not** wait on Milestone 11 for *N-light*; **11 is in scope for
G1** before product cascade.

**Scope lock:** this 5×6 is **HM450, first-65 536-locus, gene-linked**
(~2 646 genes / 51 375 CpG columns). Shared `φ`/`ρ` are gene-invariant
**architecturally**; they have not been shown to transfer to the other
~17k graph genes, EPIC, or ONT. Do **not** treat the checkpoint as the
~20k-gene product. Plan:
[`milestone-12b-full-gene-panel.md`](milestone-12b-full-gene-panel.md).
Cascade 5×6 (native P2-G) should use the **G1** panel, not copy this prefix.

**Protocol:** 5 outer folds × up to 6 restarts; orientation-aligned scores
([ADR 0008](../adr/0008-score-identifiability.md)); no TBS
([ADR 0009](../adr/0009-drop-tbs-scores.md)).

**Split (frozen):** `hub-nine-pack-5fold-v1` (34 234 samples, 470 studies).
Do **not** reuse `hub-nine-pack-3fold-v1` for this OOF.

**Product readout for every trait:** frozen nested elastic-net on encoder
scores (`mbs_enet_nested`; cascade also `rbs_enet`). Joint `mbs_e2e` is
encoder supervision / diagnostic only.

## Arms — N-light first

1. **Light (now)** — N-light@64 one-hop, mean pool, gene-linked CpGs.
   Runner: `scripts/run_12_nlight_oof.sh` (GPU 2, `batch_size` 1024 with
   VRAM probe). Run ids `stage0-12-nlight-oof-f{fold}-r{restart}`.
   Config: **30 epochs**, early-stopping patience **15**, aux heads for
   **9 disease classes + 9 cancer types** with nine-pack disease-tissue
   **n>200** (exact 200 excluded). Nine-pack 3-fold nested already
   **0.368 / 9.88 / 0.803** vs e2e **0.308 / 14.73 / 0.852** — do not 5×6
   e2e-only.
2. **Cascade (THEN — after GATE G1–G4)** — native **P2-G** unless fair **10e**
   flips the recipe. On the **G1** gene panel (12b recipe ± M11), not another
   65k-prefix 5×6. Rank vs N-light under nested enet.

**Encoder aux (n>200, this OOF):** nine-pack `disease tissue` classes above
200 samples, with pack-matched `control` as negatives (adjacent-normal
unknown). Locked lists live in `configs/experiment/stage0_12_nlight_oof.yaml`.

- **Disease (9):** Alzheimer’s 945, schizophrenia 536, SLE 341, Parkinson’s
  333, RA 250, MS 228, psoriasis 211, stroke 204, childhood asthma 202.
- **Cancer types (9):** glioma 481, HNSCC 437, AML 437, prostate 401,
  stomach 390, breast 343, melanoma 327, HCC 241, sarcoma 239.
- **Excluded at n=200:** glioblastoma, ccRCC, cervical, CRC, LUAD,
  endometrial, TCC, KIRP, LUSC, thyroid; intellectual disability 200.

**Freeze-reuse (n≥600, not extra OOF arms):** BMI, ancestry, pack-level
cancer/disease probes, optional post-hoc Alzheimer’s probe. Do **not** train
on naive pack `cancer_mask` / `disease_mask`.

## Panel / cohort / platforms

- **This OOF:** HM450 nine-pack, **65k column prefix**, then explicit
  gene-linked CpGs (**2 646 genes**). Same encoder family as M10.
- **Architecture:** ragged Deep Set with shared `φ`/`ρ` (no gene IDs in the
  encoder). That is a **capability**, not a completed whole-genome /
  cross-platform claim.
- **Not this run:** full 482k / ~20k-gene panel (12b); mixed EPIC/ONT;
  Milestone **11** fold-selected sparse panel.

## Prerequisites

| Required for N-light 5×6 | Required for cascade 5×6 | Optional / parallel |
|--------------------------|--------------------------|---------------------|
| Milestone **9** gene-only selection | 10e 1-fold S1–S4 smoke | Milestone **11** Stage B |
| N-light@64 + 3-fold nested enet | P2-G nested enet ranking | Freeze-reuse n≥600 traits |
| 30-ep budget + n>200 aux heads | | |
| One-hop smokes (done) | | |
| `hub-nine-pack-5fold-v1` | | |

Narrative: [ADR 0007](../adr/0007-crossfit-prerequisites.md),
[`post-v0-scientific-programme.md`](post-v0-scientific-programme.md).

**Next:** finish this 65k 5×6, then
[`milestone-12b-full-gene-panel.md`](milestone-12b-full-gene-panel.md)
(full gene-linked ~20k-gene sampler; cascade OOF on that graph). Expression
continue/finetune remains Milestone **13**.

## Completeness gate (added 2026-09-08)

Neither `run_12_nlight_oof.py` nor `run_12_cascade_oof.py` validated, before
declaring the campaign done, that all `folds x 6 restarts` jobs actually have
*valid* (present AND finite) primary metrics -- an existing `metrics.json`
was treated as "completed" even if its nested-enet subprocess (run with
`check=False`) failed silently. `scripts/check_12_oof_completeness.py` is the
missing read-only gate: run it against a run-prefix before reporting a 5x6 as
finished. Also documented (not yet fixed, doesn't need to be): checkpoint
*selection* uses `validation_tissue_macro_f1_then_age_mae`, computed before
nested enet ever runs -- `primary_evaluation: mbs_enet_nested` describes the
reported readout, not what picked the checkpoint. See the NOTE comments in
both `stage0_12_nlight_oof.yaml` and `stage0_12_cascade_oof.yaml`.

**12b / GATE:** product recipe (within-gene CpG sampler + minibatch gather) is
**G1** and correctly sequenced *after* this 5×6. **G2** (CpGPT/positional),
**G3** ([`milestone-12c-platform-robustness.md`](milestone-12c-platform-robustness.md)),
and **G4** (fair 10e + 10d) also block cascade. **Probes started** on GPU 0
(dense full-width N-light; CpGPT 65k ablation) — do not treat those as recipe
acceptance. See inventory in
[`milestone-12b-full-gene-panel.md`](milestone-12b-full-gene-panel.md).
Does not block the in-flight N-light run.
