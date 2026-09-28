# Milestone numbering index (Stage 0)

**Canonical checklist:** [`../TODO_PIPELINE.md`](../TODO_PIPELINE.md).

After Milestone **7F**, lettered suffixes (`7G`, `7G′`, `7H`, …) became hard to
navigate. **Current numbering uses integers 8, 9, 10, …** with optional `a/b/c`
sub-tracks. Older plan filenames and `stage0_7g_*` / `run_7g_*` artifact IDs
are **historical aliases** — do not rename on-disk run trees.

**Execution order ≠ numeric order.** Agents follow the **NOW → GATE → THEN**
block at the top of [`../TODO_PIPELINE.md`](../TODO_PIPELINE.md) (GPU queue +
gene-holdout). Do not restart closed light benchmarks.

## Map (old → new)

| New | Status | Title | Historical alias |
|----:|--------|-------|------------------|
| **1–6** | done | Bootstrap → flat → hierarchical v0.1 | unchanged |
| **7A–7F** | done | Release → nine-pack → architecture → CV → cascade topology | unchanged |
| **7E′** | done | Hub multitask + hygiene | sub-letter retained (done) |
| **8** | done | Methylation-only full eval + tissue probe | **7G** |
| **9** | done | Gene-only architecture selection (ATS) | **7G′ Stage A** (+ screens) |
| **9a** | done | DeepRVAT Tier-1 / gene-only probe grid | 7G′ Stage A screen |
| **9b** | done | Matched 16-epoch promotion (2×2 pooling retained) | 7G′ 16-ep |
| **9c** | done | Age-primary seed-mask (not adopted) | 7G′ seed-mask |
| **9d** | done | Typed-RBS R0–R5 CPU (neural typed not promoted) | 7G′ typed-RBS |
| **10** | campaign (topology locked) | Pretrained MBS/RBS scale; residual: GPU queue + 10d | **7H** |
| **10a** | done | Nine-pack grid; P2-G locked cascade finalist | 7H Track A |
| **10a-cap** | done | N-light rho_hidden=64 adopted | 7H capacity |
| **10a-warm** | done | Vector LP-FT warm-start (both pooling recipes) | 7H warmstart |
| **10a++** | done | One-hop seed-mask + multi-seed smokes | 7H correctness |
| **10b** | done | ATS seed-43 pooling | 7H Track B |
| **10c** | **done** (pack-level) | Freeze-reuse cancer/disease/BMI probes | 7H Track C |
| **10d** | pending | Reference checkpoint(s) — after cascade OOF finalist | 7H Track D / Phase 4 |
| **10e** | **done (FAIL)** | Fair S1→S2→S3→S4 — native P2-G retained; **GATE G4** | pre-OOF recipe |
| **11** | ATS done; G1 open | Fold-selected panel / gene-holdout — **GATE G1** | **7G′ Stage B** |
| **12** | N-light **done**; cascade **blocked** | 65k-prefix 5×6 closed; product cascade waits GATE | historical **Milestone 7** |
| **GATE** | G1/G3 open; G2 YES (N-light); G4 10e done | Blocks cascade OOF | see TODO_PIPELINE |
| **12b** | **← NOW** (science) | Gene-holdout + DeepRVAT sampler (**GATE G1**) | after 12 65k OOF |
| **12c** | pending | Platform robustness (**GATE G3**) — CpG/platform-mask dropout; EPIC later | after 12 65k OOF |
| **13** | deferred | Expression aux (continue/finetune **after** OOF; download first) | **7G″** |
| **14** | deferred | Optional Stage 1+ layers (a–f only) | historical **§8 Optional** |

**G2** (CpGPT DNA-sequence-embedding probe) lives under 12b probe host + TODO gate table;
no separate milestone number. **YES for N-light**; cascade smoke still in GPU queue.

## Plan file aliases

| New plan stub | Points at (do not delete) |
|---------------|---------------------------|
| [`milestone-8-methylation-eval.md`](milestone-8-methylation-eval.md) | `milestone-7g-methylation-eval.md` + tissue probe |
| [`milestone-9-gene-only-architecture.md`](milestone-9-gene-only-architecture.md) | `milestone-7g-prime-*.md` family |
| [`milestone-10-pretrained-mbs-rbs.md`](milestone-10-pretrained-mbs-rbs.md) | `milestone-7h-pretrained-mbs-rbs-campaign.md` |
| [`milestone-10e-staged-rbs-mbs-training.md`](milestone-10e-staged-rbs-mbs-training.md) | GATE G4 staged recipe (N-light unblocked) |
| [`milestone-11-fold-selected-panel.md`](milestone-11-fold-selected-panel.md) | Stage B; G1 gene-set option |
| [`milestone-12-final-oof.md`](milestone-12-final-oof.md) | ADR 0007 + programme § final OOF |
| [`milestone-12b-full-gene-panel.md`](milestone-12b-full-gene-panel.md) | GATE G1 gene utilization |
| [`milestone-12b-gene-holdout.md`](milestone-12b-gene-holdout.md) | G1 gene-holdout plumbing (random + chrom; N-light + cascade) |
| [`milestone-12c-platform-robustness.md`](milestone-12c-platform-robustness.md) | GATE G3 platform dropout / EPIC path |
| [`milestone-13-expression-auxiliary.md`](milestone-13-expression-auxiliary.md) | `milestone-7g-double-prime-expression-auxiliary.md` (+ download) |

Legacy redirects (old numbers):
[`milestone-12-expression-auxiliary.md`](milestone-12-expression-auxiliary.md) → **13**;
[`milestone-13-final-oof.md`](milestone-13-final-oof.md) → **12**.

## Agent rule

When writing new docs or commit messages, prefer **Milestone 8 / 9 / 10…**.
When citing older reports or scripts, keep the historical path/ID and add
`(alias: Milestone N)` once.
