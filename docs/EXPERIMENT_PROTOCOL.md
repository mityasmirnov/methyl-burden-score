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

Development protocol (Milestone **7E**):

```text
3 outer folds
2 restarts per fold
```

Final Stage 0 protocol (Milestone **7**, after 7A–7E):

```text
5 outer folds
up to 6 restarts per fold
```

Do not launch final 5×6 until catalog census, nine-pack matrices, multi-path
architecture, Level-1 normalization, and 7E selection are complete
([ADR 0007](adr/0007-crossfit-prerequisites.md)). A 3-fold / 1-restart smoke of
existing machinery is plumbing only—not Milestone 7.
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
| M01 | gene MBS only | flat or hier gene | — |
| M02 | gene + direct CpG | multi-path | — |
| M03 | gene + RBS + TBS + direct | multi-path | — |
| N01 | M* + Level-1 robust z | selected topology | — |

M* / N* arms are Milestone **7E** (independently trained). Also required:
transparent gene/region mean and elastic-net baselines; parameter-matched
flat vs hierarchical (same width/activation/dropout/norm); CpGPT as a
**separate** ablation. Attention-only pooling is deferred until this matrix is
complete. Report **macro-F1**, **balanced accuracy**, RMSE, R², correlations,
AUROC/AUPRC, and calibration—not accuracy/MAE alone—for selection. Residual
eval must not use an ordered prefix of holdout samples (v0.1 used first 512).

## Negative controls

### Label permutation

Shuffle phenotype labels within appropriate study strata. Performance must return to chance.

### Static-only

Remove methylation values while retaining static vectors and graph structure. Strong performance indicates locus, coverage, vocabulary, platform, or study leakage.

### Coverage-only

Use gene presence and observed-CpG counts without beta values.

### Metadata-only

Predict from study, platform, and tissue metadata alone (no methylation). Strong
performance is a confounding ceiling, not a burden-model success.

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

For each sample, average only checkpoints from folds and restarts that excluded
its full group, after **orientation alignment** ([ADR 0008](adr/0008-score-identifiability.md)).
Store model IDs, weights, and `score_polarity` in the score manifest. This is
Milestone **7** after architecture selection in **7E**. Do not retrain v0.1.

Independently trained branch ablations (gene-only vs regulatory vs direct, etc.)
are required before claiming a branch is uninformative; eval-time masking of a
jointly trained residual slot is not sufficient.

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
(default for Hub and 5c configs). With `logging.auto_tensorboard` (default
**true** whenever TensorBoard logging is on), `mbs train flat` also spawns or
reuses a TensorBoard process for that run’s `tb/` directory and prints
`tensorboard_url` / `monitor_hint` in the train summary.

```text
$MBS_ARTIFACT_ROOT/runs/<run_id>/metrics.jsonl      # one JSON record per epoch
$MBS_ARTIFACT_ROOT/runs/<run_id>/tb/                # TensorBoard events
$MBS_ARTIFACT_ROOT/runs/<run_id>/tensorboard.json   # port / pid / SSH tunnel hint
$MBS_ARTIFACT_ROOT/checkpoints/<run_id>/best.pt
$MBS_ARTIFACT_ROOT/checkpoints/<run_id>/last.pt
```

Live dashboard (second SSH session; starts/reuses TensorBoard + Rich TUI):

```bash
source scripts/activate_data_environment.sh
RUN=stage0-flat-multitask-age-tissue-v1
uv run mbs monitor --run-id "$RUN"
# --no-tensorboard / --tb-port 6007 as needed
```

SSH tunnel from your laptop: `ssh -L <port>:localhost:<port> …` then open
`http://localhost:<port>`. Prefer the port from `tensorboard.json` or the
train/monitor printout — do not launch a second manual TensorBoard for the same
run (busy port errors usually mean one is already serving).

Optional raw checks: `tail -f …/metrics.jsonl`, `nvidia-smi`,
`ps … | rg 'mbs train flat'`.

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
