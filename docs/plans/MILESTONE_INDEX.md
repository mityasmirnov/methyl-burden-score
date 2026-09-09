# Milestone numbering index (Stage 0)

**Canonical checklist:** [`../TODO_PIPELINE.md`](../TODO_PIPELINE.md).

After Milestone **7F**, lettered suffixes (`7G`, `7G′`, `7H`, …) became hard to
navigate. **Current numbering uses integers 8, 9, 10, …** with optional `a/b/c`
sub-tracks. Older plan filenames and `stage0_7g_*` / `run_7g_*` artifact IDs
are **historical aliases** — do not rename on-disk run trees.

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
| **10** | **in_progress** | Pretrained MBS/RBS scale campaign | **7H** |
| **10a** | done | Nine-pack grid; P2-G locked cascade finalist | 7H Track A |
| **10a-cap** | done | N-light rho_hidden=64 adopted | 7H capacity |
| **10a-warm** | done | Vector LP-FT warm-start (both pooling recipes) | 7H warmstart |
| **10a++** | done | One-hop seed-mask + multi-seed smokes | 7H correctness |
| **10b** | done | ATS seed-43 pooling | 7H Track B |
| **10c** | **done** (pack-level) | Freeze-reuse cancer/disease/BMI probes | 7H Track C |
| **10d** | pending | Reference checkpoint(s) — **GATE G4** after fair 10e + OOF | 7H Track D / Phase 4 |
| **10e** | **in_progress / reopen** | Fair S1→S2→S3→S4 smoke — truncated ≠ reject; **GATE G4** | pre-OOF recipe |
| **11** | pending (G1) | Fold-selected panel — in scope for **GATE G1** gene-set; does not block N-light 65k | **7G′ Stage B** |
| **12** | **in_progress ← NOW** | N-light 5×6 GPU2 (HM450 **65k-prefix** validation) | historical **Milestone 7** |
| **GATE** | pending | Post–N-light G1–G4 **blocks cascade OOF** | see TODO_PIPELINE |
| **12b** | pending | Gene utilization (**GATE G1**) — DeepRVAT sampler; probes ≠ recipe | after 12 65k OOF |
| **12c** | pending | Platform robustness (**GATE G3**) — CpG/platform-mask dropout; EPIC later | after 12 65k OOF |
| **13** | deferred | Expression aux (continue/finetune **after** OOF; download first) | **7G″** |
| **14** | deferred | Optional Stage 1+ layers (a–f only) | historical **§8 Optional** |

**G2** (positional / CpGPT probe) lives under 12b probe host + TODO gate table;
no separate milestone number.

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
| [`milestone-12c-platform-robustness.md`](milestone-12c-platform-robustness.md) | GATE G3 platform dropout / EPIC path |
| [`milestone-13-expression-auxiliary.md`](milestone-13-expression-auxiliary.md) | `milestone-7g-double-prime-expression-auxiliary.md` (+ download) |

Legacy redirects (old numbers):
[`milestone-12-expression-auxiliary.md`](milestone-12-expression-auxiliary.md) → **13**;
[`milestone-13-final-oof.md`](milestone-13-final-oof.md) → **12**.

## Agent rule

When writing new docs or commit messages, prefer **Milestone 8 / 9 / 10…**.
When citing older reports or scripts, keep the historical path/ID and add
`(alias: Milestone N)` once.
