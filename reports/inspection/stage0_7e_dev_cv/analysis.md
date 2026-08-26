# Milestone 7E analysis report

**Question:** On the frozen Age/Tissue/Sex Hub cohort, which *architecture* should
Milestone 7 use to turn a variable set of observed CpGs into gene-level (and
optional non-gene) scores?

**Selection winner (neural bake-off rule):** `N-multipath-l1a`
(late-fusion of gene + regulatory + tile + direct features, without Level-1 z).

**Read this first.** Neural Deep Set models were trained for only **2 epochs** on
the first **8 192 CpG columns** of the matrix. Linear and boosting models below
fit to convergence on the same folds and the same locus prefix. They are therefore
*not* a fair “neural vs gradient boosting” bake-off of fully trained models.
They *are* a fair check of whether methylation M-values already carry tissue/age
signal under study-held-out evaluation.

---

## 1. Data

| Item | What it is |
|------|------------|
| Cohort name | **ATS** = Age / Tissue / Sex |
| Matrix id | `matrix-hub-age-tissue-sex-full-v1` |
| Samples | 13548 GEO microarray samples (GSM) |
| Studies | GEO series (GSE). Splits never put the same study in train and test. |
| Genome | GRCh38 |
| Platforms | Illumina Infinium arrays (mostly 450K / EPIC) mixed in one matrix |
| Annotation graph (gene arms) | `graph-grch38-gencode38-five-role-v1` |
| Annotation graph (RBS/TBS) | `graph-grch38-gencode38-cgi-tile-v2` |
| Split pack | `hub-ats-7e-3fold-v1` — 3 outer study-grouped folds × 2 neural restarts |
| Locus budget in this CV | First **8 192** CpG columns (a compute ceiling, not a scientific claim that later CpGs are useless) |

Split uses **3** outer folds over **327** studies and **13548** samples.

**We do not use** sample IDs, study IDs, or platform IDs as features in methylation
models. The metadata-only control *does* use study and platform, on purpose, to
show how much phenotype can be guessed from batch labels alone.

---

## 2. What we tried to predict

Each sample can have none, some, or all of these labels. Missing labels are
**unknown**, not “healthy” or “control”.

| Target | Type | Why it matters | Metric we trust |
|--------|------|----------------|-----------------|
| **Age** | Continuous (years) | Epigenetic clocks are the usual benchmark | **MAE in years**, also R². Ignore neural `age_rmse` ≈ 1 — that is on *standardized* age, not years. |
| **Tissue / cell type** | Many classes (blood, brain, …) | Strong methylation signal; easy to cheat via study) | **Macro-F1** and **balanced accuracy** (every class counts equally) |
| **Sex** | Two classes | Sanity check; X/Y methylation | Accuracy / AUROC when both sexes are in the holdout |

ROC curves apply to **sex** (binary) and to **tissue** only in the
one-versus-rest sense (one tissue vs all others). Age does not have an ROC
curve; that would require turning years into a yes/no disease which we did not do.

---

## 3. Glossary (every abbreviation in this report)

| Term | Plain language |
|------|----------------|
| **MBS** | Methylation Burden Score — one number (or vector) per gene, built from the CpGs that map to that gene |
| **RBS** | Regulatory Burden Score — same idea for non-gene regulatory regions (CGI / cCRE-like) |
| **TBS** | Tile Burden Score — same idea for intergenic 50-CpG tiles. We do **not** assign those CpGs to the nearest gene |
| **CpG** | A cytosine-guanine dinucleotide; the usual methylation site on Illumina arrays |
| **beta value** | Fraction methylated at a CpG, roughly 0–1 |
| **M-value** | `log2(beta / (1 − beta))`. More Gaussian; standard for linear models |
| **Deep set / DeepRVAT-style** | Neural net that pools a *variable-length* list of CpGs per gene (order does not matter) then predicts phenotypes from the gene scores |
| **Flat** | One pooling step: CpGs → gene score |
| **Hierarchical (hier)** | Two pooling steps: CpGs → region → gene, plus leftover (“residual”) CpGs |
| **Level-1 (L1)** | Fold-fitted robust z-score of M-values: `(M − median) / (1.4826 × MAD)`, fit on the **training fold only**. Channel **A** = off, **B** = on |
| **MAD** | Median absolute deviation — a robust spread estimate |
| **CpGPT** | Frozen static DNA-language embeddings of each locus (not sample-specific methylation) |
| **Late fusion** | Train branches separately, then glue their *features/scores* with a linear head. In this bake-off the fusion layer used **region means**, not saved neural MBS matrices |
| **Direct** | Per-CpG elastic-net (no gene pooling) whose predictions are concatenated as extra columns |
| **OOF** | Out-of-fold — a sample is never scored by a model that trained on it |
| **3×2** | 3 outer folds × 2 random restarts (neural only) |
| **5×6** | Planned final protocol (Milestone 7); **not** this report |
| **Macro-F1** | F1 score averaged across classes; rare tissues count as much as blood |
| **Balanced accuracy** | Mean of per-class recalls |
| **MAE** | Mean absolute error (years for age) |
| **RMSE** | Root mean squared error. **Only comparable when the unit is the same** |
| **AUROC** | Area under the ROC curve; 0.5 = coin flip, 1 = perfect ranking |
| **SVA** | Surrogate Variable Analysis — unsupervised axes that capture batch. Here: **PCA on train M-values**, then residualize |
| **HGB** | Histogram Gradient Boosting (sklearn). Same family as **LightGBM**; we did not add a LightGBM dependency |
| **Elastic-net** | Linear model with both L1 (sparse) and L2 (shrinkage) penalties |
| **Ridge** | Linear model with only L2 penalty |
| **GSM / GSE** | GEO sample / series ids |
| **ATS** | Age-tissue-sex pack |
| **Hub** | EWAS Data Hub compiled packs |
| **ADR** | Architecture decision record in `docs/adr/` |

---

## 4. Schemas — what went into each model

```
Sample
  └─ observed CpGs (beta)
        ├─ optional M-value
        ├─ optional Level-1 z (train-fold MAD)
        └─ optional CpGPT static vector (locus, not sample)

Gene path     CpGs in gene regions  → pool → gene scores → phenotype heads
RBS path      CpGs in CGI/regulatory tiles → pool → RBS scores
TBS path      CpGs in intergenic tiles     → pool → TBS scores
Direct path   CpG M/z matrix → elastic-net phenotype predictions
Late fusion   [gene means | RBS means | TBS means | direct preds] → linear heads

M-value classical (this add-on)
  beta → M (8 192 loci) → Ridge / elastic-net / HGB
                         → or PCA-SVA residual M → Ridge / logistic
```

| Arm | Inputs combined | Pooling | Head |
|-----|-----------------|---------|------|
| T-mean-gene | Presence-aware mean beta per gene | Mean | Ridge age + logistic tissue |
| T-enet | Same gene means | Mean | Elastic-net age + logistic tissue |
| N-flat-gene-* | beta, M, [z], [CpGPT] on gene edges | Max Deep Set | Joint multitask linear heads |
| N-hier-gene-* | Same + region types + residual CpGs | Hierarchical Deep Set | Joint multitask linear heads |
| N-rbs / N-tbs | Same features restricted to RBS or TBS edges | Flat Deep Set | Joint multitask linear heads |
| N-gene-direct-* | Gene means + direct elastic-net preds | Late linear | Ridge / logistic |
| N-multipath-* | Gene + RBS + TBS means + direct preds | Late linear | Ridge / logistic |
| C-mvalue-ridge | M-values, 8 192 loci | None (CpG matrix) | Ridge + logistic |
| C-mvalue-enet | M-values, 8 192 loci | None | Elastic-net + logistic |
| C-mvalue-hgb | M-values, 8 192 loci | Trees | Histogram gradient boosting |
| C-mvalue-sva | M-values residualized on 10 PCA SVs | None | Ridge + logistic |
| C-metadata | Study id + platform (no methylation) | — | Ridge + logistic |

Neural encoder (flat and hierarchical, matched): GELU, dropout 0.1, LayerNorm,
CpG hidden size 64.

---

## 5. Protocol (how leakage was blocked)

1. Freeze **one** 3-fold study-grouped split (`hub-ats-7e-3fold-v1`). Every arm
   reuses it. No sample, donor, or study is in both train and the held-out test
   of the same fold.
2. Level-1 medians / MADs, PCA surrogate variables, scalers, and linear heads
   are fit on the **training studies of that fold only**.
3. Winner rule for neural architecture arms: highest mean tissue **macro-F1**,
   ties broken by age error in **years** (MAE). Transparent and metadata models
   are ceilings, not candidates for Milestone 7’s *score architecture*.

---

## 6. Results

| Model | Family | Tissue macro-F1 | Tissue balanced acc. | Age MAE (years) | Age R² | Cells |
|---|---|---:|---:|---:|---:|---:|
| Metadata-only (study + platform) | control | 0.659 ± 0.093 | 0.660 | 9.76 ± 0.92 | 0.716 | 3 |
| Late-fusion gene+RBS+TBS+direct, no Level-1 | neural | 0.329 ± 0.048 | 0.366 | 11.49 ± 1.11 | 0.624 | 6 |
| Late-fusion gene+RBS+TBS+direct, Level-1 | neural | 0.329 ± 0.048 | 0.366 | 11.49 ± 1.11 | 0.624 | 6 |
| M-value SGD elastic-net / logistic | classical | 0.324 ± 0.027 | 0.351 | — | — | 3 |
| Gene-mean elastic-net (transparent) | transparent | 0.322 ± 0.033 | 0.351 | 20.22 ± 1.94 | 0.011 | 3 |
| M-value PCA-SVA + ridge / logistic | classical | 0.290 ± 0.056 | 0.308 | 19.03 ± 4.05 | -0.503 | 3 |
| M-value ridge age + SGD-L2 logistic | classical | 0.283 ± 0.032 | 0.310 | 10.77 ± 2.85 | 0.635 | 3 |
| Late-fusion gene + direct, no Level-1 | neural | 0.270 ± 0.016 | 0.309 | 14.01 ± 1.49 | 0.462 | 6 |
| Late-fusion gene + direct, Level-1 | neural | 0.270 ± 0.016 | 0.309 | 14.01 ± 1.49 | 0.462 | 6 |
| Gene-mean linear (transparent) | transparent | 0.236 ± 0.053 | 0.276 | 15.59 ± 3.99 | 0.363 | 3 |
| Hierarchical neural, gene, Level-1, no CpGPT | neural | 0.225 ± 0.056 | 0.262 | 18.86 ± 2.52 | 0.058 | 6 |
| RBS neural branch, Level-1 | neural | 0.218 ± 0.015 | 0.284 | 18.10 ± 1.25 | -0.125 | 6 |
| Flat neural, gene, Level-1, no CpGPT | neural | 0.169 ± 0.034 | 0.227 | 16.71 ± 1.71 | -0.048 | 6 |
| Hierarchical neural, gene, Level-1 | neural | 0.160 ± 0.049 | 0.186 | 19.68 ± 2.02 | 0.009 | 6 |
| Flat neural, gene, Level-1 | neural | 0.138 ± 0.058 | 0.196 | 18.38 ± 1.46 | -0.199 | 6 |
| M-value histogram gradient boosting | classical | 0.103 ± 0.061 | 0.139 | 11.13 ± 1.34 | 0.647 | 3 |
| RBS neural branch, no Level-1 | neural | 0.102 ± 0.052 | 0.148 | 19.88 ± 0.88 | -0.397 | 6 |
| Hierarchical neural, gene, no Level-1 | neural | 0.079 ± 0.031 | 0.103 | 22.61 ± 3.18 | -0.320 | 6 |
| Flat neural, gene, no Level-1 | neural | 0.055 ± 0.028 | 0.093 | 20.99 ± 1.33 | -0.560 | 6 |
| TBS neural branch, Level-1 | neural | 0.049 ± 0.013 | 0.090 | 18.79 ± 1.81 | -0.341 | 6 |
| TBS neural branch, no Level-1 | neural | 0.020 ± 0.016 | 0.043 | 20.00 ± 2.04 | -0.489 | 6 |

Figures:

- Model schema (what is combined): `figures/model_schema.png`
- Tissue F1 and age MAE bars: `figures/arm_bars.png`
- Tissue F1 by fold: `figures/fold_heatmap.png`
- Tissue one-vs-rest ROC (M-value HGB, fold 0): `figures/roc_tissue_ovr.png`
- Sex ROC (M-value HGB, fold 0): `figures/roc_sex.png`

Classical M-value models used **8192** loci (matrix has **482379** columns). M-values on the same 8192-locus prefix as the neural bake-off. This is all ATS samples, not a sample subset. The Hub matrix has more CpG columns; using every column would be tens of GB and is not required for an architecture-matched comparison. Penalised linear models use SGD with L2 or elastic-net penalties (coordinate-descent ElasticNet / SAGA logistic did not finish on 8192 loci).

### How to read the table

- **Metadata-only** is a *confounding ceiling*: if methylation cannot beat it,
  the model may be picking up study identity rather than biology. It is expected
  to look strong because many GEO series are single-tissue, single-age-band.
- **Neural RMSE ≈ 1** for flat/hier/RBS/TBS is **not** ~1 year. Those RMSE
  values are on standardized age. Compare age using **MAE (years)** only.
- **N-multipath-l1a and l1b are nearly identical** because late fusion used the
  same region-mean features; Level-1 did not enter that fusion matrix.
- **C-mvalue-enet age is blank on purpose.** SGD elastic-net for years exploded
  (unbounded linear predictions). Tissue logistic from the same family is kept;
  do not read a trillion-year MAE as a scientific result.

---

## 7. What was missing or interrupted (honest gaps)

The 90-cell bake-off **did finish** (3 folds × arms). These gaps are about
*evaluation quality*, not a crashed trainer:

1. **Under-trained neural nets.** `max_epochs: 2` and `max_loci: 8192` were a
   compute ceiling. A fully trained flat/hier model could close the gap to
   M-value boosting. Do **not** conclude “trees beat Deep Sets” from this table.
2. **Late fusion is not neural MBS fusion.** Independent gene/RBS/TBS nets were
   trained, but the reported multipath numbers are linear models on
   **presence-aware region means** (+ direct elastic-net predictions). Saving
   per-sample score matrices and fusing *those* is still outstanding.
3. **T-mean-region** (region-mean transparent arm) was in the plan and not a
   separate named cell (gene-mean covers the transparent story).
4. **No LightGBM package.** Histogram gradient boosting is the same algorithm
   family (leaf-wise histogram trees). Installing LightGBM would not change the
   qualitative conclusion at 2-epoch neural budget.
5. **Not every CpG column.** The Hub matrix has more loci than 8 192. Using the
   full column set is a larger job (memory). The classical comparison is
   matched to the neural prefix, which is the honest architecture test.
6. **Neural ROC.** Stored neural `auroc` fields are not a 47-class tissue ROC;
   they come from a binary helper inside the training loop. Trust the HGB ROC
   figures for ranking plots.
7. **SVA is PCA-SVA.** Full Bioconductor `sva` (moderated t, iteratively
   estimated surrogate count) was not run. Ten train-only principal components
   removed as covariates is the usual first-order substitute.
8. **Sex** is incomplete in the neural summary dump (heads were trained;
   the merged table focused on tissue + age).

---

## 8. Recommendation

- **For Milestone 7’s score *topology*:** keep **multi-path** (gene + RBS + TBS
  + a direct CpG branch) with **late fusion**. That is the only architecture
  that uses noncoding tiles without stuffing them into the nearest gene
  (ADR 0006).
- **For the phenotype head:** a linear or boosted head on concatenated branch
  scores is currently stronger than 2-epoch joint DeepRVAT heads. Revisit after
  a longer neural train, still on these **same frozen folds**.
- **Do not** treat metadata-only as a model to ship. It is the leakage alarm.
- **Do not** start 5×6 OOF (Milestone 7) until 7E′ Hub disease/cancer heads
  exist and unknown labels stay unknown.

---

## 9. Files

| Path | Content |
|------|---------|
| `summary.json` / `summary.md` | Raw 90-cell bake-off dump |
| `arm_means.json` | Fold-averaged table used here |
| `classical_baselines.json` | M-value ridge / enet / HGB / PCA-SVA |
| `figures/` | Bars, heatmap, schema, ROC |
| `milestone-7e-dev-cv.canvas.tsx` | Interactive Cursor canvas source (open beside chat from the IDE canvases folder; this file is the git-tracked mirror) |
| `configs/experiment/stage0_7e_bakeoff.yaml` | Arm matrix |
| `artifacts/splits/hub-ats-7e-3fold-v1/` | Frozen folds (not in git) |
