# Stage 0 experiment protocol

## Primary questions

1. Can a shared Deep Set scorer learn useful sample–gene methylation scores from variable observed CpG sets?
2. Does typed CpG-to-region-to-gene hierarchy improve over direct CpG-to-gene aggregation?
3. Do static CpGPT sequence-adapter features improve generalization beyond raw methylation and structured annotations?
4. Are improvements preserved when entire studies and platforms are held out?
5. Are scores stable under technical replication and simulated manifest downsampling?

## Pilot data

Start with a small, auditable CpGCorpus release:

```text
2–3 blood age studies
2–3 multi-tissue studies
1–2 EPIC studies
1 technical-replicate study
1 complete study reserved as a final holdout
```

Prefer studies with QCDPB processing and clear metadata. Keep study-provided processed matrices out of the first scientific comparison unless QCDPB coverage is insufficient.

## Split policy

Grouping priority:

```text
donor and technical replicate
> cohort
> study or GSE
> sample
```

No group may cross outer folds. Inner validation is also study-grouped.

Development protocol:

```text
3 outer folds
2 restarts per fold
```

Final Stage 0 protocol:

```text
5 outer folds
up to 6 restarts per fold
```

## Training tasks

### Age

- regression on age standardized with train-fold statistics;
- report MAE in years, RMSE, Pearson correlation, Spearman correlation, and R²;
- stratify by study, tissue, platform, and age range.

### Tissue

- multiclass cross-entropy;
- report macro-F1, balanced accuracy, per-class precision/recall, and confusion matrix;
- exclude classes represented by only one study from claims of cross-study generalization.

## Baseline matrix

| ID | Representation | Topology | Pooling |
|---|---|---|---|
| B00 | gene mean beta | direct gene | mean |
| B01 | region mean beta | typed regions | mean/mean |
| B02 | max robust deviation | direct gene | max |
| B03 | MethylSPW-style linear capsule | typed regions | weighted sum |
| D00 | beta only | flat Deep Set | max |
| D01 | beta + raw NTv2 | flat Deep Set | max |
| D02 | beta + CpGPT adapter | flat Deep Set | max |
| H00 | beta only | hierarchical | max/max |
| H01 | beta + annotations | hierarchical | max/max |
| H02 | beta + CpGPT adapter + annotations | hierarchical | max/max |
| H03 | beta + MethylGPT token prior | hierarchical | max/max |
| H04 | fused static priors | hierarchical | max/max |
| H05 | CpGPT adapter + annotations | hierarchical | gated/max |
| H06 | CpGPT adapter + annotations | hierarchical | mean/mean |

Attention-only pooling is deferred until this matrix is complete.

## Negative controls

### Label permutation

Shuffle phenotype labels within appropriate study strata. Performance must return to chance.

### Static-only

Remove methylation values while retaining static vectors and graph structure. Strong performance indicates locus, coverage, vocabulary, platform, or study leakage.

### Coverage-only

Use gene presence and observed-CpG counts without beta values.

### Random locus vectors

Replace MethylGPT token priors with fixed random vectors of equal dimension.

### Shuffled graph

Permute locus-to-gene edges while preserving region size, gene size, and platform-coverage distributions.

### Random 49K panel

Compare the MethylGPT-selected panel with coverage-matched random loci.

## Robustness tests

### Manifest downsampling

For each test sample, score:

```text
all available probes
450K overlap
EPIC overlap
random 50%
random 25%
random 10%
```

Report:

- per-gene MBS correlation;
- median absolute score difference;
- phenotype-prediction change;
- fraction of genes becoming unobserved.

### Technical replicates

Report intraclass correlation and absolute agreement for MBS and phenotype predictions.

### Unseen genes

Exclude a gene subset from all phenotype heads during training and apply the shared scorer afterward. This tests transfer to new genes.

### Feature-source dropout

For fused CpGPT and MethylGPT features, randomly remove each static branch during training to prevent total dependence on one source.

## Model selection

Select checkpoints by mean normalized validation loss across tasks. Do not select on the final holdout study.

Record:

- optimization steps;
- effective observed-CpG count;
- samples by task;
- gradient norms;
- task metrics;
- study and platform metrics;
- missing-gene fraction.

## Out-of-fold score generation

For each sample, average only checkpoints from folds and restarts that excluded its full group. Store model IDs and weights in the score manifest.

## Minimum success criteria

Stage 0 does not require a specific biological benchmark value before implementation. It requires:

- a functioning exact flat baseline;
- a functioning hierarchical model;
- leakage-free study-held-out evaluation;
- reproducible static-feature export;
- stable score artifacts;
- all negative controls behaving as expected;
- at least one hierarchical configuration that is competitive with transparent baselines.

## Reporting

Every experiment report includes:

```text
Git commit
resolved configuration
canonical data release
annotation graph release
static-feature artifact IDs
fold-assignment hash
checkpoint hashes
hardware and software summary
aggregate and stratified metrics
negative-control results
known study overlap with foundation-model pretraining
```
