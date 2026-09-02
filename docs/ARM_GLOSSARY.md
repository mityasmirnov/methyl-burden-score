# Evaluation arm glossary (7E–7G′)

Canonical source: `src/mbs/inspection/arm_glossary.py`. Inspection reports embed a subset at the top of each `analysis.md`.

## Naming prefixes

| Prefix | Meaning | Examples |
|--------|---------|----------|
| **P** | **Probe** ablation (7G tissue investigation or 7G′ `-G` gene-only) | `P2-end2end-tissue-weight`, `P4-G` |
| **C-mvalue-** | **Classical** sklearn/HGB on M-values (no neural encoder) | `C-mvalue-enet`, `C-mvalue-enet-G` |
| **N-** | **Neural** encoder (CascadeDeepSet or flat/hier variants) | `N-cascade-l1`, `N-light-type` |
| **T-** | **Transparent** fixed-feature linear baselines | `T-mean-region`, `T-mean-gene` |
| **-G** | **Gene-linked** CpG panel only (`gene_cols`; 7G′ Stage A) | `P2-G`, `C-mvalue-ridge-G` |
| **-S** | Fold-safe **Selected** sparse panel (7G′ Stage B) | `C-mvalue-enetS`, `N-cascade-S` |

## Evaluation modes (cascade metrics)

- **`mbs_e2e`**: End-to-end **`MultitaskHeads`** on MBS only (no orphan/direct columns) — **Stage A primary metric**.
- **`mbs_linear_probe`**: CPU linear probe fit on **saved MBS matrix** only (representation check, not fusion).
- **`fusion_full`**: Late fusion on **`[orphan_rbs | mbs | direct_contrib]`** columns.
- **`fusion_mbs_direct`**: Late fusion on **`[mbs | direct_contrib]`** — orphan RBS ablation.

## Quick reference (most common arms)

| Arm | What it is |
|-----|------------|
| `P0-baseline` | 7G product cascade replay; late fusion on 65k loci |
| `P2-end2end-tissue-weight` | Cascade retrain with heavy tissue loss; **late fusion** test |
| `P2-G` | Gene-only cascade (P2 weights); primary metric **`mbs_e2e`** |
| `P4-G` | Gene-only cascade with mean/mean pooling; **`mbs_e2e`** |
| `P5-G-max` / `P5-G-mean` | Gene-only cascade, 30-epoch cap + early stop |
| `P2-orphan-ablation` | Full assignment; compares fusion with vs without orphan RBS |
| `C-mvalue-enet` | Elastic-net on all 65k prefix loci (7G classical tissue winner) |
| `C-mvalue-enet-G` | Same enet on **gene-linked** panel (Stage A fair comparator) |
| `C-mvalue-enetS` | **Stage B:** fold-selected sparse panel → enet |
| `N-cascade-l1` | 7F/7G product topology + late fusion |
| `N-cascade-S` | **Stage B:** locked cascade on fold-selected panel |
| `N-light-type` | FlatDeepSetRegion (lighter encoder) |
| `N-full` / `N-mbs-direct-only` | Stage B full fusion vs orphan ablation |
| `T-mean-region` | Region-mean M-values → linear heads (strong transparent baseline) |

For the full machine-readable table, run:

```bash
uv run python -c "from mbs.inspection.arm_glossary import render_full_glossary_doc; print(render_full_glossary_doc())"
```
