# Pre-OOF training recipe: staged RBS → MBS + frozen-feature heads

> **Status:** `pending` interpretation lock (2026-09-08). **Do not launch
> Milestone 12 5×6 OOF on joint `mbs_e2e` as-is.** Cheap 1-fold recipe smoke
> first.
>
> Parent: [`milestone-10-pretrained-mbs-rbs.md`](milestone-10-pretrained-mbs-rbs.md).
> Evidence: [`../../reports/inspection/stage0_7h_nine_pack_smoke/vector_vs_scalar.md`](../../reports/inspection/stage0_7h_nine_pack_smoke/vector_vs_scalar.md),
> ATS table in [`../TODO_PIPELINE.md`](../TODO_PIPELINE.md),
> [`../../reports/inspection/stage0_7g_gene_only_probe/analysis.md`](../../reports/inspection/stage0_7g_gene_only_probe/analysis.md).

## Why joint e2e is the wrong OOF recipe

Three independent findings, all on honest outer test:

1. **Region RBS (frozen + classical) beats gene MBS e2e.** Nine-pack
   `rbs_linear_probe` on ~15k region scores: best tissue **0.368**, best age
   **9.64**. Same encoder’s `mbs_e2e`: 0.333–0.355 tissue, 13.4–15.4 age.
   `mbs_linear_probe` (classical on *pooled* MBS) still beats e2e age
   (P2-G 11.57 vs 13.43). So the encoder’s **region features are good**; the
   jointly trained multitask head + max-pool gene hop **throw signal away**.
2. **Freeze / LP-FT is how we recovered vector.** Cold vector lost; warm-start
   from a frozen scalar encoder closed most of the gap. Training a random
   `gene_rho` jointly with the encoder was the failure mode — same class of
   mistake as joint e2e heads.
3. **Elastic-net >> e2e, especially N-light.** ATS `N-light-gene-mean`:
   `mbs_e2e` age MAE **17.1** vs `mbs_enet_nested` **10.3** (tissue 0.378 vs
   0.387). Nested enet is the product readout for the light model. Nine-pack
   N-light@64 e2e (0.308 / 14.7) has **not** yet been scored with enet —
   that probe is a prerequisite, not an afterthought.

**Interpretation.** OOF should train **encoders** the way DeepRVAT trains
burden models (stage the pooling hops, freeze, then heads), and **report**
frozen-feature elastic-net (`rbs_enet` / `mbs_enet` / nested) as co-primary
with e2e. Shipping 5×6 of the current joint-e2e recipe would lock in the
weaker readout.

## Staged encoder recipe (cascade)

Cheap-validate on **1 fold, ~8–15 ep**, nine-pack, age/tissue/sex only.
Promote to 3-fold only if stage-3 test beats P2-G joint e2e **and** matches
or beats that run’s own `rbs_linear_probe` / `mbs_enet`.

| Stage | What is trained | Frozen | Pooling / hop | Why |
|------:|-----------------|--------|----------------|-----|
| **S1** | CpG→region encoder (`phi`, `region_rho`) | — | **mean** CpG→region and region→(pre-MBS) so every region gets a dense gradient | Max-pool starves `region_rho` (current P2-G limitation). This is the in-flight dense stage-1 queue. Prefer **vector region embeddings** here (capacity), not scalar RBS yet. |
| **S2** | Scalar-ize / burden-ize RBS: max (or learned) over regions **inside** the gene for a scalar RBS, still **no gene MBS hop** — or keep vector RBS | S1 encoder frozen 2–4 ep, then low-LR | max is the *product* aggregation (“is there an impaired region?”); do not start with max | Matches LP-FT: new pooling module must not wreck S1. |
| **S3** | Gene MBS: `gene_rho` / learned pooler over (vector or scalar) RBS → MBS. **Not** a hard max-only hop if a small MLP/attention beats it on the 1-fold smoke | RBS stack frozen | gene hop | MBS is a compression of RBS; train it on frozen good RBS. |
| **S4** | Unfreeze all, fine-tune @ ~3e-4 | — | keep S2/S3 pooling | Same recipe that made warm vector competitive. |

**S1 vector → S2 scalar** is the default hypothesis (dense vector RBS, then
scalar burden for deployment). **S3 on vector RBS** (region_hidden MBS) is
the alternate if S2 scalar loses the RBS probe. Decide from the 1-fold
smoke, not from the cold 5-combo.

Dense stage-1 (`run_7h_dense_stage1_queue.sh`) is **S1 only** plus a naive
transplant into max/max — necessary but **not** the full S1–S4 recipe.

## Readout policy (cascade and N-light)

| Readout | Role before / during OOF |
|---------|--------------------------|
| `rbs_linear_probe` / `rbs_enet` | **Co-primary** for cascade (region features). Already computed. |
| `mbs_enet` / `mbs_enet_nested` | **Co-primary** for MBS product, **required** for N-light (ATS nested enet is the dramatic age win). Turn `include_mbs_enet` **on** for nine-pack refs that skipped it. |
| `mbs_e2e` | Secondary until S4 heads catch enet; keep for on-device / no-sklearn deploy. Improving e2e (deeper heads, separate age head, not more joint CE+Huber soup) is a **head** problem, not a reason to skip enet. |
| `C-mvalue-enet-G` | External ceiling (ATS tissue 0.388). OOF may omit 5×6 classical but should still cite this bar. |

Do **not** pick the OOF architecture with `mbs_e2e` alone. The campaign already
showed that metric ranks P2-G joint e2e first while RBS probes and N-light
enet sit higher.

## N-light (one-hop)

No region hop — S1–S2 collapse to: train `FlatDeepSetRegion` (rho=64, mean
pool) on nine-pack, **always** score `mbs_enet_nested` (and linear probe).
e2e age is known-broken (ATS 17 y vs nested enet 10 y). OOF light arm =
encoder + **enet nested heads**, not 5×6 of e2e-only.

One-hop seed-mask / multi-seed smokes stay as correctness gates on the
*encoder*; they do not justify skipping enet.

## Data: 173k Hub vs 34k nine-pack

Nine-pack (**34 234**, HM450, split `hub-nine-pack-3fold-v1`) is the **honest
encoder-training cohort**. Catalog Hub is **~173k** samples / 1 763 studies —
use it as **frozen-score reuse**, not as a leaked joint retrain.

| Use | Cohort | How |
|-----|--------|-----|
| Train S1–S4 encoders | Nine-pack (or a documented HM450 union with the same study-grouped split) | Joint age/tissue/sex (and later masked disease/cancer if labels are honest) |
| ATS heads | ATS ⊂ nine-pack / ATS split | Frozen RBS/MBS → enet / logistic |
| BMI | 2 070 labeled in nine-pack table | Freeze-reuse regression; GPU only after recipe smoke |
| Ancestry | 1 380; no class ≥1k | Freeze-reuse CE; expect weak per-class; no OOF 5×6 until collapse ADR |
| Cancer types | cancer pack case/control ~9.1k + diagnosis strings | `label_status` case/control **and** multilabel `phenotype_value`; **not** pack `cancer_mask` |
| Individual diseases (AD, PD, stroke, …) | Hub diagnosis strings | **Census before GPU** (nine-pack Hub cases: AD **945**, schizophrenia **536**, SLE **341**, PD **333**, UC **258**, MS **228**, RA **225**, psoriasis **211**, stroke **204**, … — 28 labels). Freeze-reuse logistic vs matched `label_status=control`; skip if n≪1k |
| Blood / brain packs | catalogue `sample_type=control` | **No trait head** until pack semantics are fixed |

Encoder stays **gene-invariant** (DeepRVAT). New traits = new probes on frozen
RBS and/or MBS, optionally a tiny unfrozen head. That is how 173k and extra
packs pay off without poisoning the 34k split.

## Gate before Milestone 12

1. Finish one-hop correctness smokes (encoder only).
2. Nine-pack **enet probes** on existing P2-G + N-light@64 checkpoints
   (`include_mbs_enet` / nested) — CPU, no new GPU encoder.
3. **1-fold S1–S4 smoke** vs P2-G joint e2e and vs that smoke’s RBS/MBS enet.
4. Disease/cancer (and AD/PD/stroke if n allows) **freeze-reuse probes** on
   the winning frozen scores — confirms the 173k/pack story without 5×6.
5. Only then: Milestone **12** 5×6 on **the staged recipe + enet co-primary**,
   not on the old joint-e2e 15-ep P2-G/N-light pair alone.

Dense stage-1 (GPU 2) informs S1; it does **not** skip steps 2–4.
