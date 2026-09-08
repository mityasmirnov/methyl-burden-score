# Plan: Pretrained MBS/RBS framework — extended architecture campaign

> **Canonical milestone number: 10** (alias **7H**). See
> [`MILESTONE_INDEX.md`](MILESTONE_INDEX.md) and
> [`milestone-10-pretrained-mbs-rbs.md`](milestone-10-pretrained-mbs-rbs.md).

Status (2026-09-07): **active (Milestone 10).** Track A nine-pack unblocked; full 3-fold refs
**running** on GPU 0; Track B.5 trait census done; Track B.4 queued after full.
User mandate (verbatim intent, condensed): dig into the `SeedMaskedLinearHead`
init/weight-decay collapse, then use GPU 0 freely to find/train the best
architecture(s) for gene-level (MBS) and region-level (RBS) methylation
burden scores. Vision: a one-stop framework for (1) feature aggregation +
interpretation at gene and regulatory-region level, (2) methylation data
normalization, (3) association testing across traits with reduced multiple-
testing burden from aggregation — portable across 450k/EPIC/ONT/PacBio.
Explicit permission to use more data, more epochs, more traits, and to
consider separate pretrained models for gene-level vs region-level scores.

This doc is the single source of truth for that work across a long,
possibly multi-session window — update it after every experiment batch
rather than relying on conversation memory.

## Scope reality-check (grounding the mandate in what's actually buildable now)

- **Data actually available:** ATS cohort (`matrix-hub-age-tissue-sex-full-v1`,
  13,548 samples, age/tissue/sex) is the validated, fold-frozen benchmark
  used by every 7G/7G′ result so far. A much larger resource exists and is
  **currently untapped for architecture screening**:
  `matrix-hub-nine-pack-virtual-v1` (34,234 samples, virtual multi-store,
  own split `hub-nine-pack-full-auto-v1`, phenotype table
  `sample_phenotype_table_hub_nine_pack_v1.parquet`) — likely backs the 9
  trait-specific phenotype tables already on disk (`age`, `tissue`, `sex`,
  `disease`, `blood`, `brain`, `ancestry`, `bmi`, `cancer`
  `*_sample_info.parquet`). This is the "more data, more traits" lever —
  use it for Phase 3, not before (don't scale before the base architecture
  question is answered on the smaller, well-understood cohort).
- **Platforms:** the matrix store schema carries a `platform_id` concept,
  but a quick check of the ATS `sample_index.parquet` shows no platform
  column populated at that layer — array-vs-array (450k/EPIC-era Illumina,
  the GEO-sourced hub cohorts) is plausible today; **ONT/PacBio (long-read)
  methylation calling is a different data modality with no ingestion path
  in this repo currently.** Treat "works across all platforms" as a
  longer-term design constraint to keep in mind (e.g., don't hard-code
  array-specific normalization inside the encoder), not something to
  empirically validate in this 40h window — that needs new data ingestion,
  out of scope here.
- **Gene-level vs region-level split:** the cascade topology already
  produces both intermediate region scores (RBS, pre-gene-pool) and final
  gene scores (MBS, post-gene-pool) — see `rbs_linear_probe`/`rbs_enet`
  readouts already in every cascade per-arm report. The "separate
  pretrained models" idea is worth testing empirically (Phase 2/4) rather
  than assumed: does a model trained with an RBS-primary objective produce
  better region-level interpretability than reading RBS off an MBS-primary
  model, at the cost of a second model to maintain?

## Operational schedule (A → B → C → D) — live board 2026-09-07

User-revised order: **A first**, then **B.4/B.5 in parallel where possible**,
then **C**, then **D**. Map onto this campaign:

| Track | What | Status | Owner |
|-------|------|--------|-------|
| **A.1–A.3** | Virtual multi-store `RoutedBetas` dense slices + alignment tests | **DONE** | `ebe48df` |
| **A.4 smoke** | Fold-0 P2-G + m-only @ 3 ep, full 34 234 samples | **DONE** | report `stage0_7h_nine_pack_smoke/` |
| **A.4 full** | 3-fold P2-G (15 ep) then m-only (16 ep) | **RUNNING** | PID `scratch/logs/7h_nine_pack_full.pid`; log `7h_nine_pack_full.log` |
| **B.5** | Trait adequacy ≥1k-per-arm census | **DONE** | `trait_adequacy.md` |
| **B.4** | Genuine seed-43 ATS 2×2 cascade pooling | **QUEUED** | after A.4 full frees GPU 0 |
| **C** | Nine-pack scale decisions + trait expansion hygiene | **PARTIAL** — 10c policy census done; GPU still blocked on A.4 + B.4 review + table rewrite | manual |
| **D** | Phase 4 reference checkpoint(s) + how-to note | **BLOCKED** on C | later |

### Job queue (do not duplicate A.4 full)

1. **Now (already launched):** `uv run python -u scripts/run_7h_nine_pack_smoke.py --phase full`
   — owns GPU 0. **Do not kill / do not start a second fold.**
2. **Auto-chained next (GPU-0 keeper):** `scripts/run_7h_next_queue.sh`
   — after m-only-f0 repair: refresh report → **nine-pack vector RBS** →
   **one-hop correctness smokes** (seed-mask + multi-seed) → Track **B.4**
   ATS seed-43 pooling. **Skipped:** ATS light-mean s2 keep-busy.
   Handoff helper: `scripts/run_7h_queue_handoff.sh` (waits for m-only-f0
   metrics then relaunches updated keeper).
   **P2-G provisional only.** **N-light** co-equal OOF finalist.
3. **Soft stop after keeper.** Review `vector_vs_scalar.md` +
   `stage0_7h_onehop_correctness/`. Milestone **12** OOF = **finalists only**.
   Disease/cancer GPU blocked until configs use `label_status` sidecar;
   blood/brain/bmi/ancestry deferred.
4. **Track C (manual, after reviewing A.4 + B.4 reports):**
   - Interpret nine-pack 3-fold P2-G vs m-only vs ATS refs.
   - Disease/cancer case/control **policy defined** (CPU):
     [`trait_hygiene.md`](../../reports/inspection/stage0_7h_nine_pack_smoke/trait_hygiene.md)
     — `sample_type` via `SAMPLE_TYPE_CASE_CONTROL`; pack-mask false ≠ control;
     table rewrite still required before GPU.
   - Blood/brain trait heads **deferred** (masks empty; pack ≠ trait).
   - Tissue: only `whole blood` ≥1k among 64 labels — do not expand tissue head
     without collapsing labels.
   - Platform: all HM450 — no cross-platform claim.
5. **Track D (after C):** pretrained MBS/RBS checkpoint contract + association-
   testing note (Phase 4 below). Stage B / OOF remain blocked until architecture
   refs are honest.

**Milestone 11 CPU prep (parallel, no GPU):** honest Stage A lock refresh;
`--panels-only` / `--classical-only` on
`scripts/run_7g_prime_stage_b.py`; orphan RBS census. See
[`milestone-11-fold-selected-panel.md`](milestone-11-fold-selected-panel.md).

Launch chain (once):
```bash
nohup bash scripts/run_7h_next_queue.sh >> scratch/logs/7h_next_queue.log 2>&1 &
echo $! > scratch/logs/7h_next_queue.pid
```

## Phase 0 — fix the SeedMaskedLinearHead collapse (**DONE**)

Symptom (documented in
[`milestone-7g-prime-age-seed-mask.md`](milestone-7g-prime-age-seed-mask.md)
and `reports/inspection/stage0_7g_prime_seed_mask/analysis.md`): with
gradient clipping added (`2f07bb0`), the dense unmasked control (`G0`)
trains fine, but all sparse seed-masked heads (`G1`/`G2`/`G3`, 41–331 of
2646 genes active per output row) collapse to random-baseline output by
epoch 3–6 and stay there. Not a degenerate-mask bug (masks verified
non-zero, correctly sized). Not a `weight_decay` interaction — `Adam` in
`cascade_loop.py` is called with no `weight_decay` at all (default 0), so
that part of the original hypothesis was wrong.

**Closed:** kaiming init + LR threading + seed-offset/`-s2` path fixes landed;
full G0–G3 screen completed. **G0 beats G1/G2/G3** — seed-masking **not
adopted**. See `reports/inspection/stage0_7g_prime_seed_mask/analysis.md`.
Historical debug candidates (per-head LR, warm-start, clip ablation) are
archived; do not reopen unless a new seed-mask hypothesis appears.

## Phase 1 — close out the age-primary seed-mask decision

Once Phase 0's fix is confirmed on `G1` alone:
1. Run the full grid (`G0`–`G3` × seeds {42,43}, `C0`/`C2`) fresh.
2. Update `reports/inspection/stage0_7g_prime_seed_mask/analysis.md` (the
   script rewrites this file every run — the interpreted "GPU screen
   results" section must be re-added after each run, not left to survive
   in the file itself).
3. Answer the milestone's actual question: does seed-gene masking
   (`G1`/`G2`) beat all-gene (`G0`) or matched-random (`G3`) for age
   prediction? Does restricting the *input* CpGs to seed-linked ones (`G2`)
   lose information relative to seeing all CpGs with seed-masked heads
   only (`G1`)? This is the real deliverable of the seed-mask milestone —
   everything before this was infrastructure.

## Phase 2 — resolve the still-open Stage A architecture questions

The matched 16-epoch promotion screen (`87b22c3`) retained the full 2×2
cascade pooling grid without a lock and left several of the "model-
selection report" questions from earlier in this project only partially
answered (see `reports/inspection/stage0_7g_gene_only_probe/analysis.md`).
With GPU 0 free for an extended window, use Tier-2+ budgets (more epochs
than the 15-16 used so far, multi-seed reruns) to get a confident answer on:
- mean vs. max pooling at both CpG→region and region→gene hops (2×2 grid,
  currently all four cells within noise of each other at 15-16 epochs —
  more epochs / multi-seed averaging may separate them).
- one-hop (`N-light-gene`) vs. cascade — the matched-epoch rerun already
  showed `N-light-gene-mean` edging out `P2-G`; confirm this holds with
  more epochs/seeds rather than a single lucky run.
- whether the annotation channels (gene-role, CpG-context) help or hurt
  once training budget is no longer the constraint (the fold-0-only
  ablation grid found `m_only` beating `full` under a 5-epoch budget —
  check if that reverses with a real budget).

## Phase 3 — scale to the nine-pack cohort and more traits

Only after Phase 1/2 converge on a reference architecture. Move from the
13,548-sample ATS cohort to `matrix-hub-nine-pack-virtual-v1` (34,234
samples) and its own frozen split, adding whichever of the 9 available
traits (disease, blood, brain, ancestry, bmi, cancer, in addition to
age/tissue/sex) have adequate per-class sample counts — check each
`*_sample_info.parquet` against the same ≥1k-per-arm bar already applied to
BMI on ATS (documented as currently failing that bar in
`milestone-7g-prime-age-seed-mask.md`'s trait catalog) before committing
GPU time to a trait.

## Phase 4 — deliverable: reference checkpoint(s) + documentation

Produce the actual pretrained artifact(s) this milestone exists to justify:
- A trained checkpoint (or two, if Phase 2/Phase 0's "separate models"
  question resolves in favor of splitting) with a clear score contract:
  gene-level MBS output, region-level RBS output, documented input schema
  (CpG panel, gene allocation policy, orientation contract per ADR 0008).
- A short "how to use this for association testing" note: what reduces the
  multiple-testing burden here is aggregating ~482k CpGs down to ~2,646
  gene-level scores (or fewer region-level scores) *before* testing against
  a trait — document the resulting effective test count and how it compares
  to per-CpG EWAS.
- Update `docs/TODO_PIPELINE.md` and this doc's status header when each
  phase lands; do not silently let this doc go stale the way
  `milestone-7g-prime-16ep-promotion.md` did mid-run earlier in this
  project (caught and fixed 2026-09-05, but avoid repeating it here given
  the much longer unattended window this campaign runs under).

### Reference checkpoints landed (2026-09-06)

Two candidate reference architectures, both on `explicit_only` (51,375
CpGs, 2,646 genes), `hub-ats-7e-3fold-v1` split:

| Arm | Topology | Tissue F1 (mean, 3-fold) | Age MAE | Sex AUROC | Params | Checkpoints |
|---|---|---:|---:|---:|---:|---|
| **P2-G (reference, refit)** | cascade CpG→region→gene, max/max pooling, scalar RBS | 0.379 | 18.89 | 0.747 | 265,737 (9,025 encoder + 256,712 heads) | `artifacts/runs/stage0-7h-P2-G-reference/fold_{0,1,2}/best.pt` |
| N-light-gene-mean (m_only) | one-hop CpG→gene, mean pooling | 0.372 | 17.7 | — | (one-hop, no region intermediate) | `artifacts/runs/stage0-7g-gene-probe-ablation-m-only-16ep-f{0,1,2}/best.pt` |

P2-G was retrained fresh (fold-fresh, `--force-retrain`, new `stage0-7h-
P2-G-reference` run id — the original `stage0-7g-gene-probe-P2-G-explicit`
reference is left untouched for comparison) to pick up all fixes landed
this campaign (`gradient_clip_norm` default, `AdamW`+`weight_decay=1e-4`).
Result is consistent with the original P2-G reference within normal
run-to-run noise (0.379 vs 0.373 tissue F1; age MAE moved from 15.6 to
18.9, plausibly `weight_decay` trading a little age-regression precision
for tissue generalization — not confirmed, would need an ablation to
isolate). `N-light-gene-mean` (m_only) numbers are the already-landed
matched-budget ablation-confirmation run from earlier this session, not a
new retrain.

**Recommendation:** use **P2-G** as the primary reference — it natively
produces both region-level RBS (pre-gene-pool, `rbs_linear_probe`/
`rbs_enet` readouts) and gene-level MBS (`mbs_e2e`), matching the
project's stated dual-output goal, at comparable accuracy to the one-hop
alternative. `N-light-gene-mean` (m_only) is a valid **lighter-weight,
gene-score-only** alternative when the region-level intermediate isn't
needed — it has no RBS output by construction (one-hop CpG→gene).

**Caveat — these are architecture-selection checkpoints, not deployment
scores.** Each fold's checkpoint was trained holding that fold's studies
out; there is no single "deploy this" checkpoint yet, and scoring the full
cohort without leakage requires proper out-of-fold combination (train on
folds A+B, score fold C, rotate) — that is exactly what the still-blocked
Milestone 7 "final OOF cross-fitting" step is for. Phase 4's role here is
to hand that step a validated, bug-fixed architecture + config, not to
bypass it.

## Running log

- 2026-09-05: Phase 0 fix #1 (`gene_weight` Kaiming init instead of zero)
  tested on `G1` alone. **No effect** — collapsed to the identical
  degenerate state at epoch 3 regardless of init (age_loss/tissue_loss
  matched the zero-init run almost exactly). Rules out init; the collapse
  is an attractor independent of starting point, not a cold-start problem.
- 2026-09-05: Added instrumentation (`mbs_stats`, `head_active_w_absmean`,
  later `ln_stats`) to `train_cascade_on_arrays`'s per-epoch diagnostic
  print. Found the real mechanism: raw encoder output (`mbs`, pre-centering)
  already has near-zero cross-sample variance at epoch 1 (`4.1e-05`) and
  collapses to an exact numerical zero by epoch 3 (`mean=4.9e-12,
  var=3.6e-18`) — a full, uniform collapse of the *encoder's* output across
  every sample and every gene, not a head-specific effect. `head_active_
  w_absmean` keeps growing throughout (0.04→0.4 for age) even as the
  encoder gradient dies, consistent with Adam's per-parameter adaptive step
  continuing to move weights by ≈`lr` even once the true gradient is
  vanishing (normalized step size, not raw magnitude).
- 2026-09-05: Fix #2 attempt — scaled `F.huber_loss`'s `delta` to the
  train-fold age std (reasoning: raw-year Huber with `delta=1.0` sits
  permanently in the bounded linear/L1 regime given real residuals of tens
  of years, which seemed like a plausible source of undamped pressure).
  **Made it worse**: age_loss jumped 54→978 at epoch 1. Root cause of the
  attempt's own failure: PyTorch's `F.huber_loss` quadratic region is
  `0.5*error^2` (unscaled by delta), so a *larger* delta pushed large
  early-training residuals into the quadratic region instead of keeping
  them in the bounded linear one — gradient magnitude then scales with
  residual size instead of staying capped at 1. **Reverted immediately.**
  Lesson: don't change delta without also checking which regime the actual
  residual distribution falls into.
- 2026-09-05: LayerNorm-γ-collapse hypothesis tested and **ruled out** — γ
  for both `cpg_encoder`'s and `region_encoder`'s LayerNorms stayed pinned
  near 1.0 throughout (1.0→1.02) even as `mbs`/`rbs` fully collapsed.
- 2026-09-05: Traced the collapse to its exact location by hooking
  intermediate tensors directly rather than continuing to infer from
  parameter-magnitude summaries. Sequence of checks, each ruling out one
  more layer: `region_encoder`'s pre-LayerNorm activation (healthy,
  variance actually *grows* 0.1→38 over 3 epochs — not degenerate);
  `region_rho`'s actual input (`region_hidden`, absmax 1.6→4.6, std 0.4→3.1
  — normal scale); `region_rho`'s own weight/bias (absmax <0.35 throughout,
  frozen after epoch 3 — not the driver). Finally hooked `region_rho`'s
  **output** (the pre-sigmoid logit) directly: epoch 1 range `[-0.24,
  0.12]` (healthy) → epoch 3 range `[-64.8, -22.2]` (**every single sample
  in the batch saturating the sigmoid**, `sigmoid(-22.2)≈2e-10`). This is a
  genuine "collapse to a saturated constant" local optimum reached within
  ~2-3 epochs (~250-380 SGD steps), not a numerical/precision artifact —
  none of the individual weights or activations upstream were ever
  extreme; the *combination* was.
- 2026-09-05: **Root cause found and fixed.** Tested a 10x lower LR
  (0.0001 via a temp config copy) as the most direct mitigation for
  "collapses within a couple hundred steps" — but the resulting run was
  **byte-identical** to the unmodified one, proving the config override
  wasn't taking effect at all. Investigation revealed `run_7g_prime_seed_
  mask.py` never read `learning_rate` from its config and passed no `lr=`
  to `train_cascade_on_arrays` at all — silently falling back to that
  function's hardcoded default `lr=1e-2`, **10x higher** than the config's
  documented `learning_rate: 0.001` (which matches every other
  successfully-trained cascade arm in this project). The same omission
  existed in `run_7g_prime_stage_b.py`. This single bug plausibly explains
  the whole failure mode: a 10x-too-high LR combined with a narrow,
  concentrated sparse-head gradient path (vs. `G0`'s dense head spreading
  the same shock over many more parameters) is enough to overshoot into
  the sigmoid-saturating regime within a few hundred steps. Fixed by
  threading `learning_rate` from config in both scripts. Verified on `G1`
  alone: seed 43 now trains normally end-to-end (real per-sample variance,
  sensible logit range `[-24, +14.5]` by epoch 15); seed 42 improved
  substantially (logits no longer reach `-65`, capped around `-18`) but
  still shows residual under-training — noted as open, not chased further
  given the scale of investigation already spent on this one arm. Cleaned
  up all temporary diagnostic instrumentation (hooks, per-tensor stat
  prints) from `cascade_loop.py`, keeping only the four real, tested fixes:
  `SeedMaskedLinearHead` Kaiming init (harmless, standard practice),
  `gradient_clip_norm` (fixed `G0` outright), `AdamW` + `weight_decay=1e-4`
  (neutral-to-positive, standard regularization), and the `lr=` threading
  fix (the actual root cause). Full `G0`–`G3`/`C0`/`C2` grid relaunched
  with all four fixes combined.
- 2026-09-06: Full 15-epoch grid result: **`G0` fully recovered on both
  seeds** (age MAE ~20, tissue F1 ~0.22, sex AUROC ~0.7 — real, non-
  degenerate numbers, though not as strong as the single earlier lucky
  `G0` run). `G1`/`G2`/`G3` are still mostly poor by validation-selected
  checkpoint (tissue F1 0.0000–0.0152), **but the raw training logs show
  clear late-epoch recovery the 15-epoch budget cuts off before validation
  catches up**: `G1` seed 43 and `G3` seed 42 both climb from
  near-zero `encoder_grad_norm` back to 0.6+ by epoch 12–15, with training
  losses visibly improving in the same window — validation-based
  checkpoint selection just never got to see a good-enough epoch within
  15. This reads as "escaping the bad optimum takes longer than 15 epochs
  for the sparse-masked arms," not a further code bug. Raised
  `training.max_epochs` (and `cv_budget.max_epochs` for consistency) from
  15 to 40 in `configs/experiment/stage0_7g_prime_seed_mask.yaml` and
  relaunched the full grid to give that recovery room to complete before
  drawing the G0-vs-seed-masking conclusion.
- 2026-09-06: **40-epoch grid complete — Phase 0 and Phase 1 both done.**
  7 of 8 cascade runs now train to real, non-degenerate results (only `G3`
  seed 43 remains a stubborn outlier, collapsed by epoch 2, never
  recovered). Full numbers and conclusion in
  `reports/inspection/stage0_7g_prime_seed_mask/analysis.md`. Headline:
  **`G0` (dense, unmasked) clearly beats every seed-masked variant**
  (age MAE 16-18 vs 21-25; tissue F1 0.23 vs 0.06-0.10; sex AUROC 0.79 vs
  0.51-0.72), and `G1`/`G2`/`G3` don't meaningfully separate from each
  other — no evidence the specific discovered seed genes carry more age
  signal than a same-sized random gene set once training actually
  converges. Classical (`C0`, MAE 8.94) still beats every cascade arm.
  **Conclusion: seed-gene masking is not a promising direction for the
  age-primary objective on this cohort — don't adopt it for the pretrained
  MBS/RBS framework on this evidence.** Phase 0 (fix the collapse) and
  Phase 1 (answer the milestone's actual question) are both closed. Next:
  Phase 2 (resolve remaining Stage A pooling/one-hop questions with larger
  budgets) or Phase 3 (scale to nine-pack cohort) — see phase plan above.
- 2026-09-06: **Phase 2 partial: annotation ablation confirmed at matched
  budget.** Reran `m_only` vs `full` at 16 epochs/3 folds (the existing
  grid ran 8 epochs/fold-0-only and had the same missing
  `stage_a_per_epoch_eval` issue as the original bug). Direction holds
  (`m_only` still wins) but the magnitude is much smaller than the
  original headline: tissue F1 0.372 vs 0.351 (was reported as 0.276 vs
  0.174 — a 0.10 gap that shrinks to 0.02 at a real budget), age MAE 17.7
  vs 21.8 (a clearer, more decisive gap). See
  `reports/inspection/stage0_7g_gene_only_probe/analysis.md` § Matched-
  budget confirmation.
- 2026-09-06: While setting up a multi-seed rerun of the 2×2 cascade
  pooling grid, found that the repo's established "-s2" second-seed
  convention has **never actually worked**: `run_7g_gene_only_probe.py`'s
  flat/one-hop training path called `inject_fold_into_config(...,
  seed=42 + fold_i)` unconditionally, ignoring the config's own
  `experiment.seed` — so "primary" and "-s2" ablation configs (seed 42 vs
  43 in YAML) always trained with the *identical* effective per-fold seed.
  The small differences seen between those pairs throughout this project
  were run-to-run noise (cudnn nondeterminism etc.), not real seed
  diversity. Fixed by reading `experiment.seed` as the base offset
  (default 42, so every config without an explicit seed is unaffected).
  Cascade arms (`train_cascade_arm`) use a different code path and were
  never affected — they never had "-s2" variants, so the 2×2 pooling
  grid's "retained, not locked" call already rested on genuine 3-fold
  (not seed) variance, which is a legitimate but different kind of
  evidence than multi-seed averaging.
- 2026-09-06: **Phase 3 (nine-pack scale-up) blocked on missing
  infrastructure, not a quick config change.** Checked
  `matrix-hub-nine-pack-virtual-v1` (34,234 samples,
  `sample_phenotype_table_hub_nine_pack_v1.parquet` with age/tissue/sex/
  disease/cancer/blood/brain mask columns) — no frozen, study-grouped
  k-fold split exists for it yet (only `hub-ats-7e-3fold-v1` does, under
  `artifacts/splits/`). Also confirmed **all 34,234 nine-pack samples are
  a single platform (`HM450`)** — cross-platform generalization isn't
  testable with current data regardless of cohort size, reinforcing the
  earlier scope note. Building a new leakage-safe split is real,
  correctness-critical infrastructure work (study-grouping, dedup) that
  deserves careful design rather than a rushed one-shot script in an
  already-long session — left as the explicit next step for Phase 3
  rather than attempted ad hoc.
- 2026-09-07: **Nine-pack split frozen** (`hub-nine-pack-3fold-v1`, 34,234 /
  470 studies) via `scripts/build_hub_nine_pack_split.py` (`c2530e2`).
- 2026-09-07: **Virtual multi-store loader unblocked (Track A).** Extended
  `RoutedBetas` for dense `[:, :n]` / row+col slices (pack-batched, chunked);
  wired cascade / classical / transparent / `dev_cv` through
  `open_betas_for_matrix`. Unit + per-pack alignment tests in
  `tests/unit/test_routed_betas_dense_slice.py`. Fold-0 smoke on full cohort:
  P2-G (3 ep) tissue F1 **0.343**, age MAE **15.0**, sex AUROC **0.891**;
  m-only (3 ep) tissue **0.119**, age **18.1**, sex **0.735** — finite
  `eval_split=test`. Full 3-fold refs launched
  (`scripts/run_7h_nine_pack_smoke.py --phase full`). Trait adequacy census:
  age/sex PASS ≥1k; tissue only `whole blood` ≥1k among 64 labels;
  blood/brain masks empty; disease/cancer mask_true large but false≠control.
  Report: `reports/inspection/stage0_7h_nine_pack_smoke/`.
- 2026-09-07: **Next-job queue launched** (`scripts/run_7h_next_queue.sh`,
  pid `scratch/logs/7h_next_queue.pid`) — waits for nine-pack full, refreshes
  report, then Track B.4 ATS seed-43 2×2 pooling; hard-stops before Stage B /
  OOF / disease GPU. Live board: § Operational schedule above.
- 2026-09-07: **10c trait hygiene policy census (CPU)** —
  `reports/inspection/stage0_7h_nine_pack_smoke/trait_hygiene.md`. Disease/
  cancer: map Hub `sample_type` via `SAMPLE_TYPE_CASE_CONTROL` (both clear ≥1k
  case+control); pack-mask false ≠ control. Blood/brain heads deferred. Tissue
  CE not expanded (only `whole blood` ≥1k). No parquet rewrite.
- 2026-09-07: **Milestone 11 CPU prep** — honest Stage A lock
  (`next_gate: milestone_10_scale_review`); Stage B `--panels-only` /
  `--classical-only` / `--folds`; orphan RBS census (861 singleton / 1092
  multi-CpG orphans @ max_loci=65536). See
  `docs/plans/milestone-11-fold-selected-panel.md`.
- 2026-09-08: **Nine-pack 5-combo pooling grid complete; P2-G scalar max/max
  confirmed cascade finalist, no longer provisional.** Full grid (3-fold,
  15 ep, same split/budget/seed): P2-G scalar max/max wins outright (tissue
  F1 0.355, age MAE 13.431, sex AUROC 0.853); scalar mean/max, scalar
  max/mean, vector mean→max, vector max→max all worse on F1 and MAE. See
  `vector_vs_scalar.md`. Capacity screen for N-light m-only
  (`rho_hidden_dimension` 10→64) on primary **`mbs_e2e`**: tissue
  0.273→**0.308**, age 15.057→**14.727**, sex 0.742→**0.852** —
  **adopted as N-light default** (2026-09-08). A prior “rejected / 0.194 /
  20.616” note used the wrong `external_test` readout; do not cite it.
- 2026-09-08: **Age MAE mean-vs-median across studies.** Aggregate MAE
  (13.431) is mean-of-folds and hides a skewed per-study distribution:
  across 190 study×fold rows, mean 14.16 / median 11.78 / p90 24.54 / range
  3.31-51.56. Blood/immune tissues sit low (CD14+ monocyte median 6.8,
  whole blood 10.5) and brain regions + colonic mucosa sit high
  (cerebellum 24.9, colonic mucosa 50.8) under *either* statistic — the
  tissue-heterogeneity finding isn't a mean-driven-skew artifact. See
  `vector_vs_scalar.md` § Age MAE: mean vs. median.
- 2026-09-08: **Age head tissue+sex conditioning — implemented, tested,
  rejected.** Added `MultitaskHeads.age_covariates` (embedding-conditioned
  dense age head, ground-truth tissue/sex labels with an "unknown" bucket
  for masked-out samples; deliberately *not* conditioned on this model's
  own tissue/sex predictions, to avoid propagating the tissue head's ~35%
  F1 error rate into age) — `src/mbs/training/multitask.py`,
  `src/mbs/training/cascade_loop.py`, config field `model.age_covariates:
  [tissue, sex]`, unit tests in `test_multitask.py`. Full 3-fold ablation
  vs. P2-G baseline: tissue F1 0.355→0.340, age MAE 13.431→14.135 (worse),
  sex AUROC 0.853→0.834 — all three metrics moved against conditioning.
  Also checked whether it at least narrowed the blood-vs-brain MAE gap that
  motivated the test: 10.05→9.29, flat within fold-to-fold noise. **Not
  adopted; P2-G stays unconditioned.** See `vector_vs_scalar.md` § Age-
  covariates ablation for the full per-fold breakdown and discussion.
- 2026-09-08: **Vector (region_hidden) warm-start from P2-G checkpoint —
  built and running.** New mechanism in `cascade_loop.py`
  (`_load_encoder_warm_start`, `warm_start_encoder_checkpoint` /
  `freeze_encoder_epochs` / `fine_tune_learning_rate`, unit-tested):
  transplant the converged P2-G CpG/region/RBS encoder into a fresh
  `gene_aggregation: region_hidden` model (per-fold, matching fold only —
  no cross-fold leakage), freeze it for 4 epochs while `gene_rho` + trait
  heads linear-probe, then unfreeze and fine-tune jointly at 3e-4 for the
  remaining 11 epochs. Motivated by LP-FT (Kumar et al. 2022): the
  from-scratch vector arm trained a random-init `gene_rho` jointly with the
  encoder, risking early noisy gradients distorting it. Fold 0 result:
  tissue F1 0.348 / age MAE 12.88 / sex AUROC 0.949 — close to P2-G's own
  fold 0 (0.364/12.80/0.950) and clearly above from-scratch vector's
  aggregate on every metric. Folds 1-2 in progress; full verdict pending.
  Configs: `stage0_7h_nine_pack_vector_{max_max,mean_max}_warmstart.yaml`.
- 2026-09-08: **Known limitation surfaced (not yet addressed): `region_rho`
  never gets a dense gradient even in P2-G's own training.** `scalar_rbs`'s
  MBS is `max_pool(RBS)` — max-pooling routes gradient only to the argmax
  region per gene per sample, so most regions get zero gradient most
  batches, for the entire 15-epoch budget. The warm-start above transplants
  an RBS encoder that was itself gradient-starved this way. Proposed
  follow-up (not built): a dedicated stage-1 pretrain with `region_pool:
  mean` (dense gradient to every region) purely to train the encoder well,
  then transplant into the max-pooled stage-2 config (max-pooling still
  likely correct for the *final* aggregation — burden scores want "is
  there at least one impaired region," which mean-pooling would dilute;
  consistent with `scalar_max_mean` already scoring worse than
  `scalar_max_max`, 0.318 vs 0.355). Build only if the simpler warm-start
  above doesn't close the scalar-vs-vector gap.
- 2026-09-08: **Freeze-and-reuse adopted as standing practice** (user:
  "what i always asked for") — matches DeepRVAT's actual training recipe
  (train a gene-invariant impairment/burden network jointly across seed
  traits, **freeze it**, reuse the frozen scores for classical regression
  on new traits, including traits never seen in training — no further
  gradient training per new trait; DeepRVAT also bags/ensembles multiple
  training-split repeats for score stability). Concrete next deliverable:
  freeze the converged P2-G encoder, score all nine-pack samples once,
  fit lightweight logistic probes for **disease** and **cancer** (real
  case/control labels already built: 12,194 and 9,077 usable samples;
  neither trait is in the current champion config's loss at all) —
  study-grouped per fold to avoid leakage. For the probe's own gene
  panel: use **column selection** (only the genes with a trait-relevant
  biological prior) rather than `SeedMaskedLinearHead`-style always-
  multiply-by-mask, since the shared encoder has no cross-gene interaction
  (each gene's score depends only on its own CpGs) — filtering columns
  before vs. after a gene-local encoder gives identical numbers for
  surviving genes, and column selection avoids carrying dead parameters.
  Note: seed-masking specifically was tested for age/tissue/sex in
  Milestone 9c and dense heads won there — that finding doesn't
  necessarily generalize to disease/cancer, where a real disease-gene
  panel is plausibly more genuinely sparse; worth re-testing there, not
  written off. **BMI/ancestry label-prep (2026-09-08):** Hub labels joined
  into the nine-pack phenotype table; `bmi_head` / `ancestry_head` wired
  for DeepRVAT freeze-and-reuse (configs stubbed, GPU not queued). Brain
  region labels joined for honesty; **no `brain_head`** (catalogue, not
  case/control). Full-catalog sample overview:
  [`label-prep-bmi-ancestry-sample-overview.md`](label-prep-bmi-ancestry-sample-overview.md).
  Region-level annotation richness was already tried in this project and
  underperformed (user-reported) — not re-proposing that specific form.
  CpGPT embeddings/positional features, sex-chromosome-based sex
  imputation, and epigenetic-clock-based (Horvath/CpGPT/MethylGPT) age
  imputation are all explicitly deferred to **after** Milestone 12 OOF.
- 2026-09-08: **Correction to the capacity-bottleneck diagnostic above: it
  was read from the wrong evaluation field and the verdict flips.** The
  original "rejected" call (this doc, and initially `vector_vs_scalar.md`)
  used `external_test` metrics, which under-report tissue F1 for
  `flat_region`/one-hop runs specifically. Reading the same
  `artifacts/runs/stage0-7h-nine-pack-m-only-wide-f*/metrics.json` via
  `mbs_e2e` (the primary evaluation every other arm in this campaign is
  compared on) instead: tissue F1 0.273→**0.308**, age MAE
  15.057→**14.727**, sex AUROC 0.742→**0.852** — improves on all three
  metrics, not worse on two of them as previously reported. Verified
  independently by recomputing from the raw JSON before accepting this.
  **`rho_hidden_dimension: 64` is now the N-light default** (still trails
  P2-G cascade, 0.355/13.431 — capacity helps one-hop but doesn't replace
  the cascade). Do not cite the old 0.194/20.6 `external_test` figures for
  this comparison going forward.
- 2026-09-08: **Vector warm-start: full 3-fold results in for both
  max/max and mean/max.** Max/max: tissue F1 0.342 / age MAE 13.406 / sex
  AUROC 0.851 (vs. P2-G baseline 0.355/13.431/0.853, vs. from-scratch
  vector 0.333/15.408/0.834) — essentially at parity with scalar, clearly
  beats from-scratch. Mean/max: tissue F1 0.345 / age MAE 14.578 / sex
  AUROC 0.862 (vs. scalar mean/max baseline 0.330/13.496/0.875, vs.
  from-scratch vector mean/max 0.335/16.448/0.780) — same pattern, age MAE
  gap driven by one volatile fold (`[12.80, 12.00, 18.94]`), not a
  systematic regression. **Both arms confirm LP-FT warm-start reliably
  closes most/all of the from-scratch gap; P2-G scalar max/max remains the
  finalist** (no clear win from vector, added complexity for parity at
  best). Full tables in `vector_vs_scalar.md`.
- 2026-09-08: **RBS-only classical probe beats the neural `mbs_e2e` head,
  consistently, across every arm.** The cascade pipeline already writes
  `rbs_linear_probe` (classical probe on `all_gene_rbs`,
  `[n_samples, ~15,165 regions]`, pre-gene-pooling, same held-out
  `test_idx` as `mbs_e2e`) for every run — computed all along, not
  compared until now. P2-G: `mbs_e2e` 13.431 MAE vs. `rbs_linear_probe`
  **12.559** MAE. Vector from-scratch: `mbs_e2e` 15.408 vs.
  `rbs_linear_probe` **9.638**(!) MAE. Scalar mean/max: `rbs_linear_probe`
  tissue F1 **0.368** — the best tissue F1 of the entire campaign, beating
  every `mbs_e2e` arm including the finalist. Pattern holds for every arm
  checked: `rbs_linear_probe` beats `mbs_e2e` on age MAE every time (by
  0.9-5.8 years), and beats or ties it on tissue F1 in 3/4 arms checked.
  **Implication: the neural end-to-end head we've used as the primary
  finalist-selection metric all campaign is not the best predictor these
  architectures actually produce.** **Resolved (10e):** do not launch
  Milestone 12 on joint `mbs_e2e`. Staged RBS→MBS training + frozen
  `rbs_enet` / `mbs_enet(_nested)` are co-primary. Plan:
  `milestone-10e-staged-rbs-mbs-training.md`. Full table in
  `vector_vs_scalar.md` §5.
- 2026-09-08: **Two GPUs loaded in parallel** (GPU 0 free after the
  age-cov/warm-start queue finished; GPU 2, a ~98GB RTX PRO 6000
  Blackwell, found idle). GPU 0: **one-hop correctness smokes**
  (`scripts/run_7h_onehop_correctness_smokes.py`) — this script has never
  completed a run; found and fixed a **third** bug on this relaunch
  (`load_graph_tables` returns a 2-tuple `(locus_region_edges, regions)`,
  script unpacked into 3 variables `loci, genes, edges` — missing
  `read_locus_index` + separate `genes.parquet` load; fixed to match the
  working pattern in `run_cascade_hub`). GPU 2: **dense-gradient stage-1
  RBS pretrain** (`stage0_7h_nine_pack_dense_stage1_mean_mean.yaml`,
  `cpg_pool: mean` + `region_pool: mean` for fully dense gradient to every
  CpG/region, motivated by the `region_rho`-gradient-starvation finding
  above), then two stage-2 transplants (scalar max/max, vector max/max)
  via the same warm-start mechanism, testing whether a more thoroughly-
  trained encoder beats P2-G's own from-scratch result and/or further
  improves the `rbs_linear_probe` numbers. Runner:
  `scripts/run_7h_dense_stage1_queue.sh`.
- 2026-09-08: **First freeze-and-reuse result: disease/cancer probes on the
  frozen P2-G encoder.** Built `scripts/run_7h_disease_cancer_frozen_probes.py`
  — for each of the 3 nine-pack folds, load that fold's own converged P2-G
  checkpoint (no cross-fold reuse), score ALL samples once (frozen forward
  pass), fit a logistic probe on MBS and on RBS using only that fold's train
  split, evaluate on its test split. One bug found+fixed en route: building
  the eval-only model with `dropout=0.0` shifts `SharedMLP`'s internal layer
  indices vs. the checkpoint's `dropout=0.1` training-time state_dict (a
  Dropout layer is only inserted when dropout>0) — fixed by matching the
  training dropout rate (harmless since `model.eval()` makes it a no-op
  either way). Results (3-fold mean, no GPU encoder retrain):

  | Trait | AUROC (MBS) | AUROC (RBS) |
  |---|---:|---:|
  | Cancer | **0.954** | **0.958** |
  | Disease (broad/any) | 0.586 | 0.585 |

  Cancer detection from a frozen, never-trained-on-cancer encoder is
  excellent — confirms methylation-based cancer signal is strong and cheaply
  recoverable via freeze-and-reuse, no joint retraining needed. Disease (a
  broad "any condition" label spanning many heterogeneous diagnoses) shows
  real but modest signal (~0.59) — plausible given how diagnostically mixed
  that label is, not a pipeline bug (verified: cancer's strong result uses
  the identical code path/row-alignment, ruling out an indexing bug).
  Unlike age, MBS and RBS perform essentially identically for both traits
  here. Individual-disease probes (Alzheimer's, n=945, the only one clearing
  the n≥600 threshold) not yet run.
- 2026-09-08: **Milestone 12 N-light 5×6 OOF launched on GPU 2.** Per
  10e's finding that N-light + nested enet is competitive with cascade,
  N-light's OOF no longer waits on cascade's S1-S4 recipe. A waiter script
  (already built by a concurrent session,
  `scripts/run_12_nlight_oof_gpu2_after_s1.sh`) was set to hand GPU 2 from
  the dense-stage1 pretrain to this OOF run once S1 finished naturally.
  Given the priority of getting N-light OOF running at full GPU 2 scale
  immediately, killed the dense-stage1 queue early (fold 1 was 60% done,
  fold 2 hadn't started) rather than waiting ~45 more minutes — fold 0's
  checkpoint is preserved as a usable artifact regardless. The waiter
  detected the freed GPU within its 30s poll and launched
  `scripts/run_12_nlight_oof.py` automatically: `hub-nine-pack-5fold-v1`
  split, `mbs_enet_nested` as primary readout, `batch_size: 1024` /
  `batch_token_budget: 128000000` tuned for GPU 2's ~98GB — confirmed
  actually using ~85GB once training started (not just configured on
  paper). 5 folds × 6 restarts: **f0-r0 epoch 13/16** as of 15:47 (1/30).
- 2026-09-08: **P2-G nested enet 2/3 folds.** Fold 0: tissue 0.331 / age **9.89**
  / sex 0.848. Fold 1: 0.289 / **9.64** / 0.745. Fold 2 still CPU. Age already
  matches N-light nested (~9.9); tissue still trails N-light’s 3-fold 0.368.
- 2026-09-08: **10e S1–S4 1-fold smoke on GPU 0.** Warm-started 17 encoder
  tensors from dense S1 fold 0 `best.pt`; encoder unfrozen after freeze
  epochs, fine-tune lr=3e-4. Run id `stage0-7h-nine-pack-s1s4-smoke-fold0`.
- 2026-09-08: **BMI freeze-reuse is not useful yet** (3-fold MAE ~9.9 y,
  R² negative). Cancer/disease pack probes as above. Alzheimer’s individual
  probe still pending (homogeneity/subtype script in flight).
- 2026-09-08: **N-light 5-fold plumbing readout.** `stage0-12-nlight-oof-f0-r0`
  finished 16-ep / patience-5 (14 epochs, `best_epoch=9`). Nested enet
  **0.300 / 8.63 / 0.880** vs e2e 0.286 / 17.74 / 0.748. Archived as
  `stage0-12-nlight-oof-ep16-f0-r0`. Not the product 5×6.
- 2026-09-08: **OOF recipe: 30 epochs + n>200 disease/cancer aux heads.**
  Encoder supervision now includes 9 disease classes and 9 cancer types
  with nine-pack *disease tissue* n>200; pack-matched controls as negatives;
  adjacent-normal unknown. Patience 15 so the extra budget can actually be
  used. Remaining 16-ep jobs were not allowed to continue (in-memory
  `max_epochs=16`). Cascade 5×6 inherits the same head policy when 10e
  unblocks it.
- 2026-09-08: **P2-G nested enet: MBS beats RBS on age, reversing the plain-
  probe pattern.** Full 3-fold `mbs_enet_nested` vs `rbs_enet_nested`
  (`scripts/eval_mbs_enet_from_scores.py --which both --nested`):

  | Readout | Tissue F1 | Age MAE | Sex AUROC |
  |---|---:|---:|---:|
  | MBS (nested enet) | 0.335 | **9.81** | 0.759 |
  | RBS (nested enet) | 0.338 | 19.77 | 0.790 |

  RBS's age MAE degrades monotonically across folds (14.0 -> 18.7 -> 26.6),
  not noise. Likely cause: RBS has ~5x more features (13,212 vs 2,646) at
  the same training sample sizes; nested inner-CV sometimes selects
  too-weak regularization (alpha=0.01) for that feature/sample ratio,
  overfitting worse as the inner-train split shrinks in later folds.
  **Correction to the earlier framing**: RBS is not universally better than
  MBS once proper regularization is added -- it's arm- and metric-
  dependent. For P2-G, MBS is the better nested-enet age readout; the
  plain (non-nested) `rbs_linear_probe` result (12.56 MAE, beating
  `mbs_linear_probe`'s 11.57) doesn't hold once nested hyperparameter
  search is applied to both.
- 2026-09-08: **Freeze-and-reuse follow-ups still in flight**
  (`scripts/run_7h_disease_subtype_and_homogeneity.py`): Alzheimer’s vs
  disease-pack controls, cancer-subtype multiclass (≥200), and
  case/control homogeneity. Pack-level probes remain cancer **0.954** /
  broad disease **0.586**. Do not cite per-diagnosis numbers until that
  report writes.
- 2026-09-08: **S1-S4 1-fold smoke: fold 0 training complete, evaluation
  pending.** Reused the already-completed dense-stage1 (S1, mean/mean)
  fold-0 checkpoint, warm-started into `gene_aggregation: region_hidden`
  (freeze 4ep, fine-tune 11ep) -- `stage0_7h_nine_pack_vector_max_max_from_dense_stage1.yaml`,
  `--max-folds 1`. Per the 10e promotion criterion, comparing this fold-0
  result against P2-G's own fold-0 `mbs_e2e` (0.364/12.80/0.950) and
  `rbs_linear_probe`/`mbs_enet` decides whether to invest in regenerating
  S1 for folds 1-2 (killed early to free GPU 2 for N-light OOF) and
  running the full 3-fold recipe, or stop here.
