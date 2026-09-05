# Plan: Pretrained MBS/RBS framework — extended architecture campaign

Status (2026-09-05): **active, autonomous, ~40h GPU-0 budget granted by user.**
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

## Phase 0 — fix the SeedMaskedLinearHead collapse (in progress)

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

Applied fix #1 (`src/mbs/models.py`, uncommitted as of this writing):
`SeedMaskedLinearHead.gene_weight` was zero-initialized, so
`d(loss)/d(mbs) = grad_output @ (gene_weight * seed_mask).T = 0` exactly at
step 0 for every masked head. Switched to `nn.init.kaiming_uniform_`
(matching `nn.Linear`'s own default) then masked, so active columns start
non-zero. **Currently testing on `G1` alone** (both seeds, `--arm G1`) to
see if this alone is sufficient — first empirical check is whether
`encoder_grad_norm` still decays to ~0 by epoch 3 the way it did before.

If init alone doesn't fix it, next candidates in priority order (test each
on `G1` alone before committing to a full-grid rerun — each full 8-run grid
costs ~2h wall-clock):
1. **Per-head learning rate / separate param group** — sparse heads may
   need a higher effective LR than the shared encoder's dense parameters to
   overcome the small active-column count; Adam's per-parameter adaptive
   rate should partially handle this already, so this is a secondary guess.
2. **Warm-start from `G0`'s trained encoder** — freeze or lightly fine-tune
   the encoder (already known to reach a working state under `G0`) and only
   train the masked head fresh, isolating whether the *encoder* or the
   *masked head* is where the pathology actually lives.
3. **Disable gradient clipping for masked-head runs specifically** to check
   whether global-norm clipping is inadvertently starving the sparse
   heads' already-small gradient contribution relative to the dense age
   loss term (clip shrinks everything proportionally; if age's raw
   gradient dominates the pre-clip norm, post-clip the sparse tissue/sex
   signal could become even more negligible in absolute terms).
4. **Instrument, don't guess** — if 1–3 don't resolve it, add per-epoch
   logging of `gene_weight.abs().sum()` and the fraction of `present=True`
   for masked gene columns specifically, to distinguish "head weight died"
   from "encoder stopped producing usable signal for the seed-gene subset"
   from "present-mask degenerated for the small CpG panel."

**Do not run the full G0–G3×2-seed grid again until a `G1`-only test shows
healthy `encoder_grad_norm` through epoch 15**, to avoid burning the ~2h
full-grid cost on a guess that didn't work.

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

## Running log

- 2026-09-05: Phase 0 fix #1 (`gene_weight` init) applied, testing on `G1`
  alone. Full grid to follow only if this checks out.
