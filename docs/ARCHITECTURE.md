# Stage 0 architecture

Post–Stage 0 modules (epimutation AE, REGENIE export, ComBat-met) are outlined in
[`STRATEGIC_PLAN.md`](STRATEGIC_PLAN.md); they are not Stage 0 prerequisites.

## Objective

Learn one scalar methylation burden score for every observed sample–gene pair while sharing the scoring function across genes and training traits.

The model must accept variable sets of CpGs and remain independent of a fixed array manifest.

**Public name:** deepMAT. **Package / CLI:** `methyl-burden-score` / `mbs`
(unchanged; see [ADR 0003](adr/0003-milestone-5b-phenotype-registry.md)).

## Core contracts

For sample `s`, CpG `c`, region `r`, and gene `g`:

```text
sample-specific methylation features x[s,c]
static locus features z[c]
annotation features a[c,r]
        │
        ▼
shared CpG encoder phi_cpg
        │
        ▼
permutation-invariant pooling within region
        │
        ▼
shared region encoder phi_region
        │
        ▼
permutation-invariant pooling within gene
        │
        ▼
shared compression network rho
        │
        ▼
MBS[s,g] in [0,1]
```

The score network is shared across CpGs, regions, genes, samples, and traits. Phenotype-specific parameters exist only in downstream heads.

## Input representation

Initial CpG input:

```text
beta value
M value
optional train-fold robust deviation
static CpGPT sequence-adapter vector
structured locus annotations
region-edge annotations
```

Excluded from the shared encoder:

```text
sample ID
study or GSE ID
donor ID
technical replicate ID
platform ID
raw probe ID
raw gene ID
phenotype label
```

## Exact flat reference model

The first neural baseline mirrors the DeepRVAT gene-impairment module:

```math
h_{s,c} = \phi(x_{s,c})
```

```math
v_{s,g} = \max_{c \in C_{s,g}} h_{s,c}
```

```math
MBS_{s,g} = \sigma(\rho(v_{s,g}))
```

Elementwise maximum is fixed and permutation invariant. This baseline isolates whether a shared set scorer works before testing the regulatory hierarchy.

## Hierarchical reference model

```math
h_{s,c} = \phi_{cpg}(x_{s,c})
```

```math
v_{s,r} = \max_{c \in C_{s,r}} h_{s,c}
```

```math
u_{s,r} = \phi_{region}([v_{s,r}, e_{type(r)}])
```

```math
q_{s,g} = \max_{r \in R_{s,g}} u_{s,r}
```

```math
MBS_{s,g} = \sigma(\rho(q_{s,g}))
```

The initial gene-region roles are:

1. promoter core;
2. promoter proximal;
3. five-prime region;
4. gene body;
5. three-prime region.

CpG-island relation and regulatory annotations are orthogonal features rather than a combinatorial explosion of region types.

## Optional gated pooling ablation

MethylSPWNet motivates learned CpG weighting, but Stage 0 does not assign a free parameter to every locus. Instead, a shared gate may generate within-region weights:

```math
\alpha_{s,c,r}
= \operatorname{softmax}_{c \in r}
  q^T \tanh(W_h h_{s,c} + W_t e_{type(r)} + W_a a_{c,r})
```

```math
\bar{h}_{s,r} = \sum_{c \in r} \alpha_{s,c,r} h_{s,c}
```

The gated representation is concatenated with max pooling. The max/max hierarchy remains the reference model.

## Missing-set semantics

A gene with no observed CpGs is not low burden.

```math
MBS_{s,g}=0.5,\qquad present_{s,g}=0
```

Phenotype heads consume centered, masked scores:

```math
\widetilde{MBS}_{s,g}
= present_{s,g}(MBS_{s,g}-0.5)
```

This makes an unobserved gene contribute zero to a linear head.

## Phenotype heads

Age regression:

```math
\hat{y}^{age}_s
= x_s^T\alpha
+ \sum_{g \in S_{age}} w_g \widetilde{MBS}_{s,g}
```

Tissue classification:

```math
\hat{\mathbf{y}}^{tissue}_s
= A x_s + W(\widetilde{\mathbf{MBS}}_s \odot m_{seed})
```

Heads may include nuisance covariates. Nuisance covariates do not enter the shared MBS encoder.

### Multitask packs (Milestone 5c)

Do not train one model per Hub ZIP. Use one shared encoder and linear heads with
per-sample task masks (age MSE/Huber, tissue CE, optional disease/cancer aux).
See [`plans/milestone-5c-multitask-shared-encoder.md`](plans/milestone-5c-multitask-shared-encoder.md)
and [`../configs/experiment/stage0_flat_multitask.yaml`](../configs/experiment/stage0_flat_multitask.yaml).

## Ragged representation

Do not materialize a tensor shaped as samples × genes × annotations × maximum CpGs.

Use flat segment tensors:

```text
cpg_features          [N_cpg, D]
cpg_locus_row         [N_cpg]
cpg_sample_index      [N_cpg]
cpg_to_region         [N_edges]
region_type            [N_regions]
region_to_gene         [N_regions]
gene_to_sample         [N_gene_instances]
gene_panel_index       [N_gene_instances]
```

Static locus features are looked up after collation and are not copied into every sample record.

## Training and scoring

Development protocol:

```text
3 study-grouped outer folds
2 random restarts per fold
study-grouped validation
```

Final Stage 0 protocol:

```text
5 study-grouped outer folds
up to 6 random restarts per fold
```

Every stored training-sample MBS value is out-of-fold. Technical replicates and repeated donor measurements remain in one fold.

## Deferred architecture

Not part of Stage 0:

- dynamic CpGPT or MethylGPT token states;
- LoRA or full foundation-model fine-tuning;
- imputed CpGs treated as observations;
- attention-only pooling as the reference operator;
- enhancer-to-gene edges without a versioned evidence policy;
- intergenic burden assignment to the nearest gene;
- episignature and epivariant production models.