# Methylation Burden Score

Shared, gene-level DNA methylation burden scores (**MBS**) from ragged CpG
sets. Public model name: **deepMAT**. Package / CLI: `methyl-burden-score` /
`mbs`.

Primary data: CNCB **EWAS Data Hub**. Authoritative progress:
[`docs/TODO_PIPELINE.md`](docs/TODO_PIPELINE.md) · numbers:
[`docs/ARCHITECTURE_BENCHMARKS.md`](docs/ARCHITECTURE_BENCHMARKS.md) ·
milestones: [`docs/plans/MILESTONE_INDEX.md`](docs/plans/MILESTONE_INDEX.md).

**Poster panel** (overview for talks / wall posters — will iterate):
[`docs/figures/deepmat-poster.png`](docs/figures/deepmat-poster.png)
([PDF](docs/figures/deepmat-poster.pdf)). Architecture diagrams:
[`docs/figures/nlight-cpgpt-architecture.png`](docs/figures/nlight-cpgpt-architecture.png) ·
[`docs/figures/p2g-cascade-architecture.png`](docs/figures/p2g-cascade-architecture.png).

---

## Status (2026-10-07)

**G1 gene-holdout PASS** (matched-n). **G2 CpGPT YES** (modest, N-light n=6).
Sparse gather + epoch sampler landed; flat-loop **resume** landed. K=8 sparse
train smoke finished 6/6 epochs (batch 111); post-train scoring was killed to
free GPU0 for another project — checkpoints kept.

**Do not** restart dense full-width converge (~90 h). **Do not** launch cascade
5×6 until product seed-bank path is exercised and G3 is addressed. GPU0 is
held for cowrd until the user lifts the hold.

Full order: [`docs/TODO_PIPELINE.md`](docs/TODO_PIPELINE.md).

---

## Leaderboard — best architectures

Triplets are **tissue macro-F1 ↑ / age MAE ↓ / sex AUROC ↑**.
**Product claims use frozen nested elastic-net** (`mbs_enet_nested`), not joint
`mbs_e2e` (architecture screen only). Prefer multi-seed nested gaps; single-seed
nested Δ under ~0.5 age MAE is noise (measured spread ≈0.19).

### Locked product / topology choices

| Role | Choice | Why |
|------|--------|-----|
| Cascade topology | **P2-G** scalar RBS, max/max, `explicit_only` | Nine-pack lock; staged S1–S4 **FAIL** |
| Light encoder | **N-light@64** + CpGPT (product default) | G2 age/sex win; width sweep null |
| Train genes (product) | Multi-trait **seed-bank** (`method: seed`) | G1 PASS → encoder transfers; cheaper than 20k-gene train |
| Score genes | ~all gene-linked (~20k) full CpGs | Deploy coverage; train may use ≤K/gene |
| Classical tissue ceiling | `C-mvalue-enet-G` **0.388** | Still best ATS tissue; neural has not beaten it cleanly |

### Nested product path (nine-pack)

| Rank | Setup | Tissue | Age | Sex | Notes |
|-----:|-------|-------:|----:|----:|-------|
| 1 | **N-light@64 + CpGPT** (fold 0, n=6) | 0.355 | **8.10** ±0.14 | **0.879** ±0.015 | **G2 YES** — modest; ranges overlap |
| 2 | N-light@64 3-fold nested | **0.368** | 9.88 | 0.803 | Light path without CpGPT |
| 3 | N-light 5×6 OOF nested (30/30) | 0.316 | 9.59 | 0.822 | Closed 65k-prefix validation |
| 4 | P2-G cascade nested 3-fold | 0.335 | 9.81 | 0.759 | Cascade nested baseline |
| — | Full-width + CpGPT (~19.6k genes, n=1) | **0.388** | 11.67 | **0.962** | Best tissue/sex on record; **undertrained**; dense converge killed |

### Topology & classical (screen / ATS)

| Claim | Cite | Verdict |
|-------|------|---------|
| Cascade topology | P2-G `mbs_e2e` **0.355 / 13.43 / 0.853** | **Locked** |
| Cascade train recipe | Fair S1→S4 vs native P2-G | Staged **FAIL** → native P2-G |
| ATS cascade e2e | P2-G **0.373** | Locked Stage A on ATS |
| Seed-gene masking (9c) | G0 ≫ G1–G3 | **Not adopted** |
| Freeze-reuse cancer | AUROC **~0.954** | Useful; BMI not yet |

Full tables + caveats:
[`docs/ARCHITECTURE_BENCHMARKS.md`](docs/ARCHITECTURE_BENCHMARKS.md).

### GATEs

| ID | Status | One-liner |
|----|--------|-----------|
| **G1** | **PASS** (matched-n) | Gene-invariant encoder transfers; chrom optional only |
| **G2** | **YES** (N-light, modest) | CpGPT on for product once arch locked |
| **G3** | not started | Platform / CpG dropout robustness |
| **G4** | 10e FAIL; 10d pending | Native P2-G; checkpoint after cascade finalist |

---

## Next steps (brief)

1. **GPU0** — leave free for cowrd until the hold is lifted.
2. **Product train** — when a card is free: seed-bank flat train
   (`method: seed` + CpGPT + `sparse_betas` + `max_cpgs_per_gene` +
   `training.resume: auto`), then score ~20k genes full-CpG. Prefer this over
   all-gene full-width train (ablation only).
3. **Optional** — `reeval_only` from
   `artifacts/checkpoints/stage0-12b-sampler-smoke-k8-f0/best.pt` if a K=8
   score report is needed (train curve already logged).
4. **Still blocked** — cascade 5×6 OOF; dense full-width converge restart.
5. **Later** — G3 platform; cascade OOF on G1 panel (native P2-G + CpGPT);
   10d reference checkpoint; paper checklist (`docs/TODO_PAPER.md`).

---

## Design principles

1. DeepRVAT-style shared scorer over CpGs → typed regions → genes → traits.
2. Ragged observed CpGs (not a fixed probe panel at deploy time).
3. Product path: **RBS → gene MBS** + leftover **direct**; **no TBS** ([ADR 0009](docs/adr/0009-drop-tbs-scores.md)).
4. Study-grouped cross-fitting for reported scores (Milestone **12**).
5. Discovery CpGs select seed *genes*; training/deploy use all observed
   gene-linked CpGs ([ADR 0012](docs/adr/0012-seed-gene-discovery-vs-deployment-input.md)).

```text
Observed CpGs
  → typed regions → RBS
       → gene-allocated RBS → pool → MBS[s, g]
  → orphan multi-CpG regions → orphan RBS
  → remaining CpGs → direct
```

Locked cascade encoder: **`P2-G`** (`CascadeDeepSet`, scalar RBS, max/max,
`explicit_only`). Light encoder: **N-light@64**. Orientation:
[ADR 0008](docs/adr/0008-score-identifiability.md).

---

## Quick start (`power-horse`)

```bash
cd /data/projects/methyl-burden-score
cp .env.example .env   # first time
source scripts/activate_data_environment.sh
uv sync --all-groups --extra training --extra analysis --frozen
uv run mbs doctor
uv run pytest tests/unit
```

All artifacts stay under `/data` (`$MBS_*` roots). See
[`docs/WORKSPACE.md`](docs/WORKSPACE.md).

```bash
make doctor && make test-fast
uv run mbs catalog refresh-release
uv run mbs train flat --overfit-fixture
# Live run: uv run mbs monitor --run-id <run-id>
```

Long jobs: set `training.resume: auto` so a kill/preemption can continue from
`last.pt` (flat loop; §3.6).

---

## Docs map

| Topic | Doc |
|-------|-----|
| Live checklist | [`docs/TODO_PIPELINE.md`](docs/TODO_PIPELINE.md) |
| Poster / talk panel | [`docs/figures/deepmat-poster.png`](docs/figures/deepmat-poster.png) |
| Benchmark ledger | [`docs/ARCHITECTURE_BENCHMARKS.md`](docs/ARCHITECTURE_BENCHMARKS.md) |
| Sparse gather (§3.5) | [`docs/plans/milestone-12b-sparse-gather-sampler.md`](docs/plans/milestone-12b-sparse-gather-sampler.md) |
| Resume (§3.6) | [`docs/plans/milestone-12b-resume-checkpoint.md`](docs/plans/milestone-12b-resume-checkpoint.md) |
| Scoring / architecture | [`docs/SCORING_PIPELINE.md`](docs/SCORING_PIPELINE.md), [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) |
| Data | [`docs/DATA_CATALOG.md`](docs/DATA_CATALOG.md), [`docs/EWAS_DATA.md`](docs/EWAS_DATA.md) |
| Programme | [`docs/STRATEGIC_PLAN.md`](docs/STRATEGIC_PLAN.md), [`docs/plans/post-v0-scientific-programme.md`](docs/plans/post-v0-scientific-programme.md) |
| Key ADRs | [0002](docs/adr/0002-ewas-datahub-primary-source.md) Hub · [0007](docs/adr/0007-crossfit-prerequisites.md) OOF · [0009](docs/adr/0009-drop-tbs-scores.md) no TBS · [0010](docs/adr/0010-gene-allocation-policy.md) `explicit_only` · [0012](docs/adr/0012-seed-gene-discovery-vs-deployment-input.md) discovery vs deploy |

Nine-pack campaign board:
[`reports/inspection/stage0_7h_nine_pack_smoke/analysis.md`](reports/inspection/stage0_7h_nine_pack_smoke/analysis.md).

## Repository policy

**Committed:** source, SQL/schemas, YAML, docs/ADRs, small fixtures, inspection
reports (not raw matrices).

**Never committed:** methylation matrices, IDATs/BAMs, checkpoints/embeddings,
secrets, `artifacts/` run dumps (except documented small reports).

## Licensing

Source-code license TBD. Hub/Atlas data, papers, and pretrained weights have
separate licenses — review before redistribution.
