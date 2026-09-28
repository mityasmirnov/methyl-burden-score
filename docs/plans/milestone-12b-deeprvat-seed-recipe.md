# DeepRVAT-aligned seed recipe (GATE G1 / product encoder)

> Status (2026-09-28): **locked design**. G1 leakage probe = **random**
> holdout (DeepRVAT-like); chromosome optional / not queued. The **product**
> encoder path after architecture lock follows DeepRVAT: train φ/ρ on a
> multi-trait seed gene bank, score genome-wide (or the complement), nested
> enet post-hoc.
> Parent: [`milestone-12b-gene-holdout.md`](milestone-12b-gene-holdout.md) ·
> ADR [`0011`](../adr/0011-seed-gene-sources.md) / [`0012`](../adr/0012-seed-gene-discovery-vs-deployment-input.md).

## Question / Approaches / Results / Verdict

| Field | Content |
|-------|---------|
| **Question** | How should deepMAT use seed genes the way DeepRVAT does — without reviving 9c head masks? |
| **Approaches** | (A) Head masks — **rejected** (9c). (B) Random gene-holdout — G1 leakage probe. (C) **Train on fold seed-gene union, score complement / all genes** — DeepRVAT-aligned product path. |
| **Results** | Pending (seed-supervised arm + multi-trait bank not yet run). |
| **Verdict** | **Adopt (C) after arch lock.** Keep CpG-first discovery (best CpG p → gene). Multi-trait seed bank is the priority product step. Ensemble of scorers only after the experiment campaign closes. |

## Locked decisions

| Choice | Decision | Why |
|--------|----------|-----|
| Seed **usage** | Train φ/ρ on seed genes (union); score complement or genome-wide | DeepRVAT core recipe; 9c masks failed |
| Head masks | **Never revive** | G0 ≫ G1–G3 on 9c |
| Discovery | **CpG-first**: per-CpG association on outer-train → gene via **best (lowest) CpG p**; optional gene-level meta-p later | Fastest / simplest; user direction 2026-09-28 |
| Expand | ADR 0012: all `explicit_only` linked CpGs of seed genes | Discovery ≠ input panel |
| Multi-trait seed bank | **After arch lock**: discover seeds per eligible trait, **union** for encoder training; trait-specific readouts stay post-hoc nested enet | DeepRVAT shared scorer improves with many phenotypes; matches “traits in training once arch locked” |
| Hybrid prior | Optional `hybrid_fold` = `external_clean` Atlas ∪ `internal_fold`, remap via ADR 0010 `explicit_only` only | Prior without Atlas symbol leakage |
| Ensemble | **Deferred** until the experiment campaign finishes | Avoid multiplying compute before recipe lock |
| G1 probe | Random → nested (cascade then N-light); chrom optional | DeepRVAT-like; not the product seed recipe |

## Recipe (product, after arch lock)

```text
outer-train samples only
  → per eligible trait: univariate CpG association (p-values)
  → gene score = min_p over that gene's linked CpGs   (best CpG)
  → take top-K and/or FDR/Bonferroni survivors per trait
  → seed_bank = UNION over traits
  → (optional) seed_bank |= external_clean Atlas genes (explicit_only remap)
  → expand each seed gene → all linked CpGs (ADR 0012)
  → train φ/ρ (+ trait aux heads) on seed genes only
  → score MBS for all genes (or heldout complement for G1-style eval)
  → nested elastic-net per trait on frozen MBS (post-hoc)
```

Config sketch:

```yaml
training:
  gene_holdout:
    enabled: true
    method: seed                 # DeepRVAT-aligned (not random/chromosome)
    seed_discovery: best_cpg_p   # CpG-first; meta_p reserved
    n_genes_per_trait: 256
    traits: [age, tissue, sex]   # expand after arch lock
    # hybrid_fold:
    #   external_clean: path/to/atlas_panel.parquet
```

## Relation to G1 gene-holdout

| Arm | Train genes | Purpose |
|-----|-------------|---------|
| `random` | Random holdout | G1 leakage probe (gate) |
| `chromosome` | Chrom holdout | Optional locality probe only |
| `seed` | Fold seed-gene **union** | DeepRVAT-faithful product path |

Nested enet on heldout/complement columns remains the transfer readout. Do **not**
interpret random-holdout numbers as the product seed recipe.

## Non-goals (this brief)

- Reviving phenotype-head seed masks.
- Porting SKAT/burden/MAF Beta(1,25) verbatim (wrong domain).
- Building the full trait-universe runners before arch lock.
- Ensemble (`n_repeats` DeepRVAT-style) before the campaign closes.

## Open questions

- Top-K vs FDR threshold for best-CpG gene selection (start with top-K = 256/trait like today's panel).
- Whether product OOF scores **all** genes or only non-seed complement for the nested readout (prefer all genes for product; complement for G1 transfer metric).
