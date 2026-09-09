# Architecture & setup benchmarks (master ledger)

**Updated:** 2026-09-09  
**Role:** One-stop list of encoder setups / architectures tested and how they
compare. Prefer this file for decisions; detailed dumps stay in the linked
reports.

| Need | Go to |
|------|--------|
| Live campaign board | [`../reports/inspection/stage0_7h_nine_pack_smoke/analysis.md`](../reports/inspection/stage0_7h_nine_pack_smoke/analysis.md) |
| Pooling grid detail | [`../reports/inspection/stage0_7h_nine_pack_smoke/vector_vs_scalar.md`](../reports/inspection/stage0_7h_nine_pack_smoke/vector_vs_scalar.md) |
| Arm name glossary | [`ARM_GLOSSARY.md`](ARM_GLOSSARY.md) |
| Model / trainer map | [`ARCHITECTURE.md`](ARCHITECTURE.md) |
| Pipeline status | [`TODO_PIPELINE.md`](TODO_PIPELINE.md) |

---

## 1. Locked decisions (do not silently reopen)

| Decision | Lock | Evidence |
|----------|------|----------|
| Cascade **topology** | **P2-G** = `CascadeDeepSet`, `scalar_rbs`, **max/max**, `explicit_only` | Nine-pack 5-combo grid |
| Cascade **training recipe** (pre-OOF) | **Native P2-G** (not staged S1–S4) | Fair 10e fold-0 FAIL |
| Light encoder | **N-light@64** (`FlatDeepSetRegion`, mean pool, `rho_hidden=64`) | Capacity screen + nested enet |
| Product readout | Frozen **`mbs_enet_nested`** (+ `rbs_enet` for cascade) | Probes ≫ joint e2e |
| Age head tissue+sex covariates | **Rejected** | All 3 metrics worse |
| Cold / warm **vector** hop | Warm helps; **does not replace P2-G** | LP-FT arms |
| Seed-mask gene panels (G1) | **Not adopted** | One-hop + cascade smokes |
| TBS / tile scores | **Dropped** | ADR 0009 |
| Gene allocation | **`explicit_only`** | ADR 0010 |
| This OOF platform / panel | HM450 nine-pack, **65k-prefix** (~2.6k genes) | Not ~20k product (12b) |

**Metric convention.** Triplets are **tissue macro-F1 ↑ / age MAE ↓ / sex AUROC ↑**.
`mbs_e2e` = jointly trained neural heads (architecture **screen** only).
Product claims use frozen nested elastic-net unless noted.

---

## 2. Cohorts & setups

| Setup ID | Cohort | Split | Panel | Role |
|----------|--------|-------|-------|------|
| **ATS** | Hub age-tissue-sex full (~13 548) | `hub-ats-7e-3fold-v1` | 65k prefix, gene-linked ~51 375 CpGs | Milestone **9** architecture screen |
| **Nine-pack** | Virtual nine packs (**34 234**), HM450 | `hub-nine-pack-3fold-v1` | same 65k prefix → ~2 646 genes | Milestone **10** scale lock |
| **OOF (in flight)** | same nine-pack | `hub-nine-pack-5fold-v1` | 65k-prefix N-light 5×6 | Milestone **12** validation |
| **Catalog freeze-reuse** | Hub ~173k | not joint retrain | frozen encoder scores | Extra traits (disease/cancer/BMI) |
| **12b (pending)** | full gene-linked HM450 | TBD | ~19.6k genes + within-gene sampler | Product cascade panel |

Encoder remains **gene-invariant** (ragged Deep Set). Platform transfer
(EPIC/ONT) is **not** claimed from these numbers.

---

## 3. Encoder families (what was tested)

| Family | Aggregation | Typical arms | Trainer |
|--------|-------------|--------------|---------|
| **Cascade** (`CascadeDeepSet`) | CpG → typed region → **RBS** → gene **MBS** | `P2-G`, vector `region_hidden`, pooling grid | `mbs train cascade` |
| **N-light** (`FlatDeepSetRegion`) | CpG (+ optional annotation) → pool by gene → MBS | `N-light@64`, one-hop smokes | `mbs train flat` / Stage A light |
| **Flat / hierarchical v0.1** | frozen baselines | `deepMAT-flat-v0.1`, hierarchical-v0.1 | do **not** retrain |
| **Classical** | sklearn on M-values (no neural encoder) | `C-mvalue-enet`, `C-mvalue-enet-G` | transparent baselines |
| **Transparent means** | fixed region/gene means → linear | `T-mean-region` | baselines |

**Cascade knobs tested**

| Knob | Values tried | Winner / reject |
|------|--------------|-----------------|
| `gene_aggregation` | `scalar_rbs` vs `region_hidden` (vector) | **scalar** topology; warm vector ≠ flip |
| `cpg_pool` / `region_pool` | max/max, mean/max, max/mean, mean/mean (S1) | Product **max/max**; mean useful for dense S1 grads only |
| Warm-start / LP-FT | freeze ~4 ep → FT @ 3e-4 | Helps vector; not finalist |
| Staged S1→S2→S3→S4 | dense S1 → scalarize → gene hop → FT | **Rejected** vs native P2-G |
| Age covariates on age head | tissue+sex embeddings | **Rejected** |
| `rho_hidden` (N-light) | 10 → 64 | **64 adopted** |

---

## 4. ATS scale (Milestone 9) — tissue F1 reference

Honest outer test, `explicit_only`, ~51 375 gene-linked CpGs.

| Question | Best arm | Tissue F1 | Notes |
|----------|----------|----------:|-------|
| Classical on gene panel | `C-mvalue-enet-G` | **0.388** (±0.018) | Still the tissue ceiling here |
| Best MBS readout (same encoder) | `P2-G` **`mbs_enet`** | 0.385 (±0.053) | Elastic-net heads |
| Best cascade joint path | `P2-G` **`mbs_e2e`** | 0.373 (±0.038) | max/max, 15 ep |
| One-hop mean (16-ep) | `N-light-gene-mean` | ~0.378 e2e | Near P2; age MAE ~17 |
| Seed-mask | G0 ≫ G1–G3 | — | **Not adopted** |

Do **not** cite pre-fix ~0.67–0.70 `mbs_e2e` (train+val+test leak).

---

## 5. Nine-pack architecture screen (`mbs_e2e`, 3-fold)

Split `hub-nine-pack-3fold-v1`, 15 epochs unless noted. **Chooses topology only.**

| Arm | Setup | Tissue F1 | Age MAE | Sex AUROC | Verdict |
|-----|-------|----------:|--------:|----------:|---------|
| **P2-G scalar max/max** | cascade, scalar_rbs, max/max | **0.355** | **13.431** | 0.853 | **LOCKED cascade topology** |
| vector mean→max **warm** | LP-FT from scalar mean/max | 0.345 | 14.578 | 0.862 | Better than cold; not finalist |
| vector max→max **warm** | LP-FT from P2-G | 0.342 | **13.406** | 0.851 | Same |
| scalar mean/max | cascade | 0.330 | 13.496 | **0.875** | Loses tissue/age |
| vector mean→max cold | cascade | 0.335 | 16.448 | 0.780 | Reject cold vector |
| vector max→max cold | cascade | 0.333 | 15.408 | 0.834 | Reject cold vector |
| scalar max/mean | cascade | 0.318 | 17.426 | 0.812 | Reject |
| **N-light@64** | flat region, mean, rho=64 | 0.308 | 14.727 | 0.852 | Light default (e2e) |
| N-light@10 | historical | 0.273 | 15.057 | 0.742 | Superseded |
| P2-G + age covariates | age head embeds tissue+sex | 0.340 | 14.135 | 0.834 | **Rejected** |

Config pointers: `configs/experiment/stage0_7h_nine_pack_p2_g.yaml`,
`stage0_7h_nine_pack_vector_*`, `stage0_7h_nine_pack_m_only_wide.yaml`.

---

## 6. Product readout: nested / frozen enet (nine-pack)

Joint `mbs_e2e` **under-ranks** good encoders. Co-primary product numbers:

| Arm | Readout | Tissue F1 | Age MAE | Sex AUROC | Notes |
|-----|---------|----------:|--------:|----------:|-------|
| **N-light@64** | `mbs_enet_nested` 3-fold | **0.368** | **9.88** | 0.803 | Light product path; beats own e2e hard |
| N-light@64 | `mbs_e2e` 3-fold | 0.308 | 14.73 | 0.852 | Diagnostic only |
| **P2-G** | `mbs_enet_nested` 3-fold | **0.335** | **9.81** | 0.759 | Cascade nested baseline |
| P2-G | `mbs_e2e` 3-fold | 0.355 | 13.43 | 0.853 | Topology screen |
| P2-G | `rbs_linear_probe` (campaign) | up to **0.368** | down to **9.64** | — | Region features strong |
| P2-G fold 0 | `mbs_enet_nested` | 0.331 | 9.89 | 0.848 | Fair-smoke comparator |

**Interpretation.** Under nested enet, N-light no longer “trails cascade.”
OOF product score must not be 5×6 of joint `mbs_e2e`.

---

## 7. Staged RBS→MBS training (Milestone 10e)

**Question:** Does S1→S2→S3→S4 beat native P2-G?

| Stage | Intended | What ran in fair smoke |
|------:|----------|------------------------|
| S1 | Dense mean grads (prefer vector) | Reused scalar mean/mean fold-0 |
| S2 | Freeze; scalarize max/max; **no gene hop** | `stage0-7h-nine-pack-staged-s2-fold0` |
| S3 | Freeze RBS; train `gene_rho` | `…-staged-s3-fold0` |
| S4 | Unfreeze FT @ 3e-4 | `…-staged-s4-fold0` |

Earlier “S1–S4” job **skipped S2** (vector LP-FT transplant) — also negative.

### Fold-0 gate vs native P2-G

| Run | Readout | Tissue F1 | Age MAE | Sex AUROC |
|-----|---------|----------:|--------:|----------:|
| S1 dense mean/mean | `mbs_e2e` | 0.348 | 14.52 | 0.875 |
| Truncated LP-FT (no S2) | `mbs_e2e` | 0.300 | 14.97 | 0.900 |
| **Fair S4** | `mbs_e2e` | **0.298** | 15.32 | 0.901 |
| Fair S4 | `mbs_enet_nested` | **0.298** | 11.03 | 0.787 |
| Fair S4 | `rbs_linear_probe` | 0.322 | 10.86 | 0.919 |
| **P2-G fold 0** | `mbs_e2e` | **0.364** | **12.80** | **0.950** |
| **P2-G fold 0** | `mbs_enet_nested` | **0.331** | **9.89** | 0.848 |

**Verdict:** Fair staged recipe **FAIL**. **No staged 3-fold.** Cascade recipe =
native P2-G. Gate JSON:
[`staged_s1s4_fold0_gate.json`](../reports/inspection/stage0_7h_nine_pack_smoke/staged_s1s4_fold0_gate.json).
Plan: [`plans/milestone-10e-staged-rbs-mbs-training.md`](plans/milestone-10e-staged-rbs-mbs-training.md).

---

## 8. One-hop correctness (N-light smokes)

| Test | Result | Adopt? |
|------|--------|--------|
| Seed-mask G0 vs G1 (fold 0) | G0 **0.335 / 17.3 / 0.75** ≫ G1 **0.146 / 21.7 / 0.56** | **No** (mask collapses) |
| Multi-seed 42/43/44 | Tissue F1 span **0.027**, distinct | Correctness **pass** |

Report: [`../reports/inspection/stage0_7h_onehop_correctness/`](../reports/inspection/stage0_7h_onehop_correctness/).

---

## 9. Freeze-and-reuse (extra traits on frozen P2-G)

Encoder **not** jointly retrained. Pack-level probes on frozen scores:

| Trait | Readout (approx.) | Useful? |
|-------|-------------------|---------|
| Pack cancer (broad) | AUROC **~0.954** (MBS) / **0.958** (RBS) | Yes |
| Alzheimer’s | AUROC **0.838** | Yes |
| Broad disease | AUROC **~0.586** | Modest |
| Cancer subtype | macro-F1 **~0.202** | Weak but distinguishable |
| BMI | MAE ~9.9, negative R² | Not useful yet |

Detail:
[`disease_cancer_frozen_probes/analysis.md`](../reports/inspection/stage0_7h_nine_pack_smoke/disease_cancer_frozen_probes/analysis.md).

Encoder **aux** BCE for OOF uses disease/cancer classes with nine-pack
disease-tissue **n>200** (Milestone 12 recipe).

---

## 10. Side-by-side: what to cite for what

| Claim type | Cite | Do not cite |
|------------|------|-------------|
| Cascade **topology** winner | Nine-pack P2-G `mbs_e2e` **0.355 / 13.43 / 0.853** | Warm vector as “better cascade” |
| Light **product** path | N-light nested **0.368 / 9.88 / 0.803** | N-light e2e alone |
| Cascade **nested** baseline | P2-G nested **0.335 / 9.81 / 0.759** | Staged S4 e2e |
| Classical ceiling (ATS) | `C-mvalue-enet-G` **0.388** | Pre-fix leak numbers |
| Training recipe | Native P2-G (10e FAIL) | Truncated LP-FT as “full staged” |
| OOF completeness | Milestone 12 N-light 5×6 (in flight) | Partial fold means as final |

---

## 11. In flight / next (not in leaderboard yet)

| Track | Status | Notes |
|-------|--------|-------|
| Milestone **12** N-light 5×6 | Running GPU 2 | 65k-prefix validation; nested enet product |
| GATE G1 / **12b** | Pending | ~20k-gene sampler; product cascade panel |
| GATE G2 CpGPT / positional | Pending | Matched smokes only |
| GATE G3 / **12c** platform | Pending | Dropout / EPIC path |
| **10d** reference checkpoint | Pending | After OOF finalist |
| Cascade 5×6 | Blocked on G1–G3 (+10d) | Native **P2-G** recipe |

---

## 12. Config & run ID cheat sheet

| Setup | Config / runner | Artifact run id (examples) |
|-------|-----------------|----------------------------|
| P2-G nine-pack | `stage0_7h_nine_pack_p2_g.yaml` | `stage0-7h-nine-pack-P2-G` |
| N-light@64 screen | `stage0_7h_nine_pack_m_only_wide.yaml` | `stage0-7h-nine-pack-m-only-wide-f{0,1,2}` |
| Vector warm max→max | `stage0_7h_nine_pack_vector_max_max_warmstart.yaml` | warmstart staging |
| Dense S1 | `stage0_7h_nine_pack_dense_stage1_mean_mean.yaml` | `…-dense-stage1-mean-mean` |
| Fair staged S2–S4 | `stage0_7h_nine_pack_staged_s{2,3,4}_*.yaml` + `scripts/run_7h_staged_s1s4_fold0.py` | `…-staged-s{2,3,4}-fold0` |
| N-light OOF | `stage0_12_nlight_oof.yaml` + `scripts/run_12_nlight_oof.sh` | `stage0-12-nlight-oof-f*-r*` |

GPU policy: **GPU 2** may be fully occupied (N-light). **GPU 0** keeps spare
VRAM (batch-capped smokes). Do not use GPU 1.
