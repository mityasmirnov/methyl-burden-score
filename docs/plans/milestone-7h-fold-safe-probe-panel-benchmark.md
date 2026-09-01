# Plan: 7H fold-safe probe-panel and association-product benchmark

Status: **pending** (required before Milestone 7 final OOF).
Parents: [7G tissue investigation](milestone-7g-cascade-tissue-investigation.md),
[ADR 0010](../adr/0010-score-export-vs-phenotype-comparator.md), and
[ADR 0009](../adr/0009-drop-tbs-scores.md).

## Decision question

Does deepMAT add value over a classical M-value elastic-net when both methods
receive exactly the same phenotype-informed CpG panel, while preserving the
association product the model is intended to create?

This milestone is a **controlled benchmark**, not permission to make the
production score vocabulary phenotype-dependent.

## Why this is useful

The existing 7G comparison is fair in sample folds and 65,536-column budget but
not maximally informative about aggregation:

- every arm receives the first 65,536 canonical columns, not a biologically or
  predictively selected panel;
- `C-mvalue-enet` performs embedded shrinkage over all columns;
- a DeepSet benefits when a selected gene is represented by enough of its CpGs
  to estimate a stable within-gene set score;
- the current 7F direct artifact is one contribution per task, not a retained
  sample×CpG association block.

A fold-safe selected-and-expanded panel makes the comparison answer the intended
question: with the same informative measurements, is learned gene/region
aggregation better than a sparse linear model?

## C-mvalue-enetS definition

`C-mvalue-enetS` means **classical M-value elastic-net on a fold-selected,
gene/region-expanded panel**:

1. beta → M-value;
2. fit median imputation and scaling on outer-train samples only;
3. select seed probes inside outer train using repeated study-grouped inner CV;
4. expand seeds through evidence-backed gene and regulatory-region annotations;
5. refit the elastic-net on the expanded panel;
6. evaluate once on the untouched outer test studies.

The suffix `S` is part of the method name and must not be used for the current
unscreened `C-mvalue-enet`.

## Primary selection method

Use elastic-net stability selection rather than a single full-training fit.

For each outer fold and trait:

- age: elastic-net regression, importance `|coefficient|`;
- sex: binary logistic elastic-net, importance `|coefficient|`;
- tissue: multinomial/one-vs-rest logistic elastic-net, importance equal to the
  L2 norm across class coefficients;
- run repeated study-grouped inner splits and tune `alpha` and `l1_ratio`
  without outer-test access;
- rank each CpG by selection frequency first, then standardized coefficient
  magnitude.

Construct one common multi-trait seed panel with at most **10,000 unique CpGs**.
Use equal trait quotas initially; unfilled quota may be reassigned. Report the
age/sex/tissue overlap and results for 2,000, 5,000, and 10,000 seed ceilings.

### Association-screen sensitivity arm

Univariate association screening is useful for interpretation but is not the
primary selector because public-study tissue labels are highly confounded with
study. If run, it must:

- use outer-train data only;
- estimate effects only in studies with within-study phenotype variation;
- meta-analyse study-level effects rather than treat all samples as IID;
- control FDR within trait;
- combine traits by normalized rank, not raw p-value.

Call this sensitivity method `S-assoc`; never mix its outer-test associations
into training-panel selection.

## Evidence-backed expansion

Starting from the seed CpGs:

1. If a seed has a versioned typed gene-region edge, select the linked gene.
2. Include all compatible observed CpGs for that gene from the active
   450K/EPIC/EPICv2 union that pass the cross-study coverage threshold.
3. If a non-gene seed lies in a well-defined multi-CpG regulatory interval,
   include the other eligible CpGs from that same `region_id`.
4. Otherwise keep the seed as a direct CpG.
5. Never use unrestricted nearest-gene allocation.
6. If expansion exceeds 65,536 loci, reduce the seed ceiling before truncating
   within a selected gene. Seed CpGs always remain included.

The manifest records seed reason, trait, stability, coefficient magnitude,
gene/region expansion source, platform coverage, and final column order.

## Exact comparison

Within each outer fold, freeze one panel manifest and give the identical locus
set and sample masks to every arm:

| Arm | Input use | Output |
|-----|-----------|--------|
| `C-mvalue-enetS` | standardized M-values, no aggregation | phenotype predictions |
| `N-deepmat-selected` | CpG→typed region→gene DeepSet + direct | MBS, direct CpGs, predictions |
| `T-mean-selected` | presence-aware gene/region means | predictions |
| `N-light-type` | `[M, one-hot region type, observed]` before mean aggregation | MBS + direct |
| `N-mbs-direct-only` | omit orphan block | MBS + direct |

Use the P4/P5-winning pooling, P2 loss weights, and locked fusion solver.
Primary endpoint is held-out-study tissue macro-F1. Also report age MAE/R², sex
AUROC, per-study metrics, calibration, runtime, peak memory, feature count, and
fold-to-fold panel stability.

## Orphan RBS rules

- One qualified genomic `region_id` equals one feature.
- Never aggregate all orphan regions together and never aggregate by type.
- Require a versioned interval and at least two eligible observed CpGs.
- Regions without an evidence-backed gene edge remain orphan; do not force them
  into MBS by nearest gene.
- Singleton/unstructured loci remain direct.
- If zero regions qualify, the orphan block is empty and
  `N-mbs-direct-only` becomes the preferred lightweight product.

## Association-product contract

Final OOF export must contain:

```text
mbs.zarr                  sample × gene
gene_present.zarr         sample × gene
orphan_rbs.zarr           sample × qualified region (possibly zero columns)
direct_cpg.zarr           sample × retained direct locus, or a lossless matrix view
direct_locus_index.parquet
sample_index.parquet
score_manifest.json
```

`direct_contrib.zarr` may remain as a phenotype diagnostic but cannot substitute
for `direct_cpg.zarr` in downstream association work.

## Leakage prevention

Every phenotype-derived operation—imputation, normalization, univariate tests,
elastic-net tuning, stability ranking, seed union, gene/region expansion
thresholds, and final refit—is fitted using outer-train data only. Panel
manifests are fold-specific and hashed. Test studies are touched once.

After method selection, a deployment panel may be frozen using all development
data, but that model is not used to claim unbiased development performance.

## Acceptance

- P4/P5 configuration is locked.
- Three outer folds × two restarts for neural selected-panel arms.
- Identical hashed panel per fold across deepMAT and `C-mvalue-enetS`.
- Direct-CpG association artifact preserves locus identity.
- Orphan policy above is implemented and tested.
- Report recommends one final OOF configuration without test-fold tuning.
