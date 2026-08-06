# Stage 0 experiment protocol

## Primary questions

1. Can a shared Deep Set scorer learn useful sample–gene methylation scores from variable observed CpG sets?
2. Does typed CpG-to-region-to-gene hierarchy improve over direct CpG-to-gene aggregation?
3. Do static CpGPT sequence-adapter features improve generalization beyond raw methylation and structured annotations?
4. Are improvements preserved when entire studies and platforms are held out?
5. Are scores stable under technical replication and simulated manifest downsampling?

## Pilot data

Primary open source is **EWAS Data Hub**
([ADR 0002](adr/0002-ewas-datahub-primary-source.md)):

```text
1–few labeling GSEs from EWAS_db/ (per-sample GSM*.txt betas)
and/or a small baseline subset under download/
2–3 blood / age-oriented packs or studies when scaling beyond the first pilot
1 complete study reserved as a final holdout
```

Prefer GMQN-normalized Hub profiles with clear sample info archives. See
[`EWAS_DATA.md`](EWAS_DATA.md) and the labeling GSE list in
[`CPGCORPUS_STAGE0.md`](CPGCORPUS_STAGE0.md) (those GSEs are present under
`EWAS_db/` even when absent from CpGCorpus S3).

CpGCorpus Arrow (requester-pays) is an optional alternate if already on disk;
do not treat it as the default Stage 0 pilot path.

### Multitask Hub packs (Milestone 5c)

Do **not** train one model per `*_methylation_v1.zip`. Train one shared burden
encoder with linear age + tissue heads (disease/cancer aux optional) and
per-sample task masks. Pack roles, folder layout, and config schema:
[`plans/milestone-5c-multitask-shared-encoder.md`](plans/milestone-5c-multitask-shared-encoder.md).

## Split policy

Grouping priority:

```text
donor and technical replicate
> cohort
> study or GSE
> sample
```

No group may cross outer folds. Inner validation is also study-grouped.

For **multi-study** runs (Milestone 5b+), study-grouped holdout is mandatory:
the same `study_id` must not appear in more than one of
`train` / `validation` / `external_test`. Single-study smoke pilots (for example
GSE35069 donor-grouped 4/2) remain valid for debugging only and must not define
the final performance story.

Use `mbs.evaluation.build_study_grouped_split` and the phenotype registry
(`configs/data/phenotype_registry.yaml`) to record split roles.

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
- report **MAE** in years, **RMSE**, Pearson correlation, Spearman correlation, and R²;
- stratify by study, tissue, platform, and age range (`mbs.evaluation.regression_metrics`,
  `metrics_by_group`).

### Tissue / multiclass cell type

- multiclass cross-entropy;
- report **macro-F1**, **balanced accuracy**, per-class precision/recall, and
  **confusion matrix**;
- exclude classes represented by only one study from claims of cross-study
  generalization (`mbs.evaluation.multiclass_metrics`).

### Binary disease / condition

- report **AUROC** and **AUPRC** (`mbs.evaluation.binary_auroc_auprc`);
- always stratify by holdout study and platform.

### Cross-study generalization

- report performance **by holdout study** and **by platform**;
- never claim generalization from a single-study donor split alone.

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

## Training monitoring

Flat / multitask runs write under `$MBS_ARTIFACT_ROOT` when
`logging.tensorboard` and `logging.jsonl` are enabled in the experiment YAML
(default for Hub and 5c configs).

```text
$MBS_ARTIFACT_ROOT/runs/<run_id>/metrics.jsonl   # one JSON record per epoch
$MBS_ARTIFACT_ROOT/runs/<run_id>/tb/             # TensorBoard events
$MBS_ARTIFACT_ROOT/checkpoints/<run_id>/best.pt
$MBS_ARTIFACT_ROOT/checkpoints/<run_id>/last.pt
```

Live 5c MVP example (`stage0-flat-multitask-age-tissue-v1`):

```bash
source scripts/activate_data_environment.sh
RUN=stage0-flat-multitask-age-tissue-v1

# Preferred: live terminal dashboard
uv run mbs monitor --run-id "$RUN" \
  --config configs/experiment/stage0_flat_multitask.yaml

ps -eo pid,etime,%cpu,%mem,cmd | rg 'mbs train flat'
tail -f "$MBS_ARTIFACT_ROOT/runs/$RUN/metrics.jsonl"
nvidia-smi
uv run tensorboard --logdir "$MBS_ARTIFACT_ROOT/runs/$RUN/tb" --bind_all --port 6006
```

Do not interpret `val_accuracy == 0` on study-holdout tissue CE as model failure
when the holdout study’s tissue class is absent from train (closed-set
multiclass). Prefer `val_loss` / `val_mae` trends and final external-test
destandardized age MAE (years). Storage layer guidance for 5c:
[`plans/milestone-5c-multitask-shared-encoder.md`](plans/milestone-5c-multitask-shared-encoder.md)
§ Storage recommendations.

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
