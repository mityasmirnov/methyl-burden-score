# Pre-OOF training recipe: staged RBS → MBS + frozen-feature heads

> **Status:** `done` — **negative** (2026-09-08). 1-fold S1–S4 smoke loses to
> native P2-G; **do not** promote staged recipe to 3-fold / cascade OOF.
> Cascade finalist remains **native P2-G scalar max/max**. **N-light 5×6**
> OOF is running on GPU 2 (30-ep + n>200 aux).
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
3. **Elastic-net >> e2e, especially N-light — and it is a readout gap, not
   an encoder gap.** `mbs_linear_probe` on the **same** MBS features also
   beats `mbs_e2e`. Enet does per-trait nested-CV α/ℓ1; the neural head gets
   uniform `weight_decay=1e-4` plus a 10× weaker age loss (`0.3` vs tissue
   `3.0`) inside a non-convex 15-ep soup. Nine-pack N-light@64 **3-fold**
   `mbs_enet_nested`: tissue **0.368**, age **9.88**, sex **0.803** vs that
   arm’s `mbs_e2e` **0.308 / 14.73 / 0.852**. Fold 0 nested age **9.57** is
   the best single age number in the campaign. Under enet, N-light no longer
   “trails cascade.” **Do not try to close this with a better `mbs_e2e`
   head before OOF.** Joint multitask training produces the encoder;
   **every reported trait** (age/tissue/sex included) uses frozen-feature
   nested enet. Re-select cascade vs vector vs N-light **under enet**, not
   `mbs_e2e` — vector’s richer features may yet win.

**Interpretation.** OOF trains **encoders** with joint multitask loss, then
**always** reports frozen nested enet per trait (age/tissue/sex included).
Shipping 5×6 of joint `mbs_e2e` as the product score would lock in the
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
transplant into max/max. Architecture-wise, `gene_aggregation: region_hidden`
already **is** S3 (`gene_rho` is a 2-layer MLP over region vectors, not a
max). The running GPU-2 queue tests S1 → transplant into that hop → LP-FT,
and **skips S2** (scalar-ize RBS with no gene hop). Informative for S1; **not**
the full recipe. After S1, GPU 2 hands off to N-light OOF (transplants
cancelled).

## Readout policy (cascade and N-light)

| Readout | Role before / during OOF |
|---------|--------------------------|
| `rbs_linear_probe` / `rbs_enet` | **Co-primary** for cascade (region features). Already computed. |
| `mbs_enet` / `mbs_enet_nested` | **Co-primary** for MBS product, **required** for N-light (ATS nested enet is the dramatic age win). Turn `include_mbs_enet` **on** for nine-pack refs that skipped it. |
| `mbs_e2e` | **Encoder-training diagnostic only.** Do not deploy or pick finalists on it. Keep the numbers in the report. |
| `C-mvalue-enet-G` | External ceiling (ATS tissue 0.388). OOF may omit 5×6 classical but should still cite this bar. |

Do **not** pick the OOF architecture with `mbs_e2e` alone. The campaign already
showed that metric ranks P2-G joint e2e first while RBS probes and N-light
enet sit higher.

## N-light (one-hop) — first Milestone 12 arm

No region hop — S1–S4 collapse to: train `FlatDeepSetRegion` (rho=64, mean
pool) as a gene-invariant encoder, **always** score `mbs_enet_nested`.
Nine-pack 3-fold nested already matches or beats P2-G joint e2e on tissue
and crushes it on age. **N-light 5×6 starts first** (GPU 2, max VRAM);
cascade 5×6 still waits on an S1–S4 smoke.

One-hop smokes (done): G0 ≫ G1 (seed-mask **not** adopted); multi-seed
42/43/44 distinct (tissue F1 span 0.027). Encoder correctness gate is
cleared.

**Confirmed at nine-pack scale (2026-09-08), not just ATS.** N-light@64
`mbs_enet_nested` (all 3 folds, `stage0-7h-nine-pack-m-only-wide-f{0,1,2}`):
tissue F1 **0.368**, age MAE **9.88**, sex AUROC 0.803 — vs. this arm's own
`mbs_e2e` 0.308 / 14.727 / 0.852. Tissue F1 and age MAE both jump
dramatically (age MAE nearly matches P2-G's own nested-enet number, 9.89);
sex AUROC is the one metric that got *worse* under enet here. This closes
most of the light-vs-cascade gap once both are read via enet instead of
e2e — a real shift in the cost/performance tradeoff worth weighing before
picking Milestone 12's light finalist. Both one-hop correctness gates
(seed-mask, multi-seed) also passed cleanly on this same run (see
`reports/inspection/stage0_7h_onehop_correctness/analysis.md`).

## Data: 34k encoder, 173k freeze-reuse, cross-platform later

Nine-pack (**34 234**, **HM450**, split `hub-nine-pack-3fold-v1` for M10;
OOF uses `hub-nine-pack-5fold-v1`) is the honest **encoder-training**
cohort. Catalog Hub is **~173k** / 1 763 studies.

**Platform constraint (must not be forgotten):** the encoder is
**gene-invariant** and consumes **ragged CpG sets**. That is how it must
work on **HM450, EPIC (more probes/gene), and ONT methylation (even more
probes/gene)** without a fixed array manifest. Frozen-feature scoring only
transfers cleanly to samples the trained encoder can encode — today that
is HM450 nine-pack. Catalog 173k is **~110k HM450 / ~57k EPIC / ~6k EPICv2**;
ONT is not in this Hub table. Do **not** promise 173k as free extra data.
Do **not** bake HM450 probe IDs into the encoder. Platform-transfer eval
is a **post-architecture** gate, not extra OOF arms.

| Use | Cohort | How |
|-----|--------|-----|
| Train encoders | Nine-pack HM450 | Age/tissue/sex **plus** n>200 disease classes and cancer types (aux BCE; pack-matched controls). Not pack masks. |
| ATS / age / tissue / sex **report** | same split, frozen scores | nested enet / logistic — not `mbs_e2e` |
| BMI | 2 070 | freeze-reuse if n≥600 (yes) |
| Ancestry | 1 380 | freeze-reuse if n≥600 (yes; no class ≥1k) |
| Cancer types | nine-pack disease-tissue **n>200** | **encoder aux heads** (9 types). Remaining types / pack C/C = freeze-reuse. `label_status` + diagnosis; **not** `cancer_mask`. |
| Individual diseases | nine-pack disease-tissue **n>200** | **encoder aux heads** (9 classes). n≥600 freeze-reuse probes remain optional (AD 945). Matched `label_status=control`. |
| Blood / brain packs | catalogue / control-only | **No trait head** |
| EPIC / ONT | later matrices | same gene-invariant weights; do not retrain a new architecture |

Encoder stays **gene-invariant** (DeepRVAT). New traits = new probes on frozen
RBS and/or MBS, optionally a tiny unfrozen head. That is how 173k and extra
packs pay off without poisoning the 34k split.

## Gate before Milestone 12 (split by arm)

**N-light 5×6 (running):** GPU 2 **f0-r0 ~epoch 20/30**, patience 15,
n>200 disease/cancer aux heads (val disease ~0.76 / cancer ~0.85). 16-ep
plumbing archived (nested **0.300 / 8.63 / 0.880**). Split
`hub-nine-pack-5fold-v1`. Runner: `scripts/run_12_nlight_oof.sh`.

**Cascade 5×6 (recipe decided):**

1. Dense S1 **fold 0 done** (queue stopped mid fold 1 to free GPU 2).
2. P2-G nested enet **3/3** (MBS **0.335 / 9.81 / 0.759**; RBS age 19.8).
3. 1-fold **S1–S4 smoke done — negative** (e2e 0.300/14.97/0.900 vs P2-G
   0.364/12.80/0.950). **No staged 3-fold.** Cascade OOF = **native P2-G**.
4. Pack-level freeze-reuse **done** (cancer strong; AD 0.838; broad disease modest).
