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
| **10a++** | running | One-hop seed-mask + multi-seed smokes | 7H correctness |
| **10b** | done | ATS seed-43 pooling | 7H Track B |
| **10c** | partial | Trait labels / freeze-reuse prep | 7H Track C |
| **10d** | pending | Reference checkpoint(s) after staged recipe + enet | 7H Track D / Phase 4 |
| **10e** | **pending (blocks 12)** | Staged RBS→MBS train + frozen enet co-primary | pre-OOF recipe |
| **11** | deferred | Fold-selected panel Stage B (**parallel; does not block 12**) | **7G′ Stage B** |
| **12** | blocked | Final study-grouped OOF (5×6) — needs **10e** recipe smoke, **not 11** | historical **Milestone 7** |
| **13** | deferred | Expression aux (continue/finetune **after** OOF; download first) | **7G″** |
| **14** | deferred | Optional Stage 1+ layers (a–f only) | historical **§8 Optional** |

## Plan file aliases

| New plan stub | Points at (do not delete) |
|---------------|---------------------------|
| [`milestone-8-methylation-eval.md`](milestone-8-methylation-eval.md) | `milestone-7g-methylation-eval.md` + tissue probe |
| [`milestone-9-gene-only-architecture.md`](milestone-9-gene-only-architecture.md) | `milestone-7g-prime-*.md` family |
| [`milestone-10-pretrained-mbs-rbs.md`](milestone-10-pretrained-mbs-rbs.md) | `milestone-7h-pretrained-mbs-rbs-campaign.md` |
| [`milestone-10e-staged-rbs-mbs-training.md`](milestone-10e-staged-rbs-mbs-training.md) | pre-OOF staged recipe (blocks 12) |
| [`milestone-11-fold-selected-panel.md`](milestone-11-fold-selected-panel.md) | Stage B sections of `milestone-7g-prime-matched-probe-lightweight.md` |
| [`milestone-12-final-oof.md`](milestone-12-final-oof.md) | ADR 0007 + programme § final OOF |
| [`milestone-13-expression-auxiliary.md`](milestone-13-expression-auxiliary.md) | `milestone-7g-double-prime-expression-auxiliary.md` (+ download) |

Legacy redirects (old numbers):
[`milestone-12-expression-auxiliary.md`](milestone-12-expression-auxiliary.md) → **13**;
[`milestone-13-final-oof.md`](milestone-13-final-oof.md) → **12**.

## Agent rule

When writing new docs or commit messages, prefer **Milestone 8 / 9 / 10…**.
When citing older reports or scripts, keep the historical path/ID and add
`(alias: Milestone N)` once.
