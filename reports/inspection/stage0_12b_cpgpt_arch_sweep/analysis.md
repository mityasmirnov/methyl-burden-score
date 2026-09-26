# CpGPT architecture sweep — verdict: no capacity effect detected

Updated: `2026-09-26` — **complete, 6/6 variants** (v2, 15.0 h GPU0 wall).

## Hypothesis under test (and its outcome)

With CpGPT2M enabled the per-CpG input width jumps **24 → 152** dims while
`phi`/`rho` and `cpg_hidden` stayed at **64**. The sweep asked whether the
locked N-light shape was therefore *capacity-starved*, understating CpGPT.

**Outcome: not supported.** Across a 4× range of `phi`/`rho` (64→256), a 2×
range of `cpg_hidden` (64→128), and a depth/dropout arm, **no configuration
beat the in-sweep control by the pre-registered margin.**

## Results (ranked by nested age MAE, the primary metric)

| variant | cpg_hidden | φ/ρ | layers | dropout | nested age | nested tissue | nested sex | h |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `c64-w128` | 64 | 128 | 2 | 0.1 | **7.685** | 0.367 | 0.827 | 2.4 |
| `c128-w256-d3` | 128 | 256 | 3 | 0.2 | 7.948 | **0.372** | 0.890 | 2.7 |
| `c64-w64` *(control)* | 64 | 64 | 2 | 0.1 | 8.117 | 0.366 | 0.865 | 2.3 |
| `c128-w256` | 128 | 256 | 2 | 0.1 | 8.238 | 0.365 | 0.842 | 2.6 |
| `c128-w128` | 128 | 128 | 2 | 0.1 | 8.383 | 0.370 | 0.844 | 2.4 |
| `c64-w256` | 64 | 256 | 2 | 0.1 | 8.478 | 0.361 | **0.904** | 2.7 |

All arms: batch 128 (matched, explicit), fold 0, 65 536-locus prefix, 16
epochs, CpGPT2M on, single seed.

## The decisive comparison

| sample | n | age MAE range | width | sd |
|---|---:|---|---:|---:|
| **6 different architectures** (this sweep) | 6 | 7.685–8.478 | 0.793 | **0.267** |
| **1 fixed architecture, repeated** (multirestart + v1 control) | 7 | 7.955–8.482 | 0.527 | **0.187** |

Varying architecture across a 4× capacity range produces spread of the *same
order* as simply re-running one identical config. The sweep is, to within its
resolution, **re-measuring noise**.

Supporting detail: nested **tissue F1 varies by only 0.011** (0.361–0.372)
across the entire grid — essentially flat. If capacity were binding anywhere,
tissue is where the earlier CpGPT regression appeared, and it does not move.

## Against the pre-registered bar

The bar (set *before* results, from measured fixed-config spread): a single-seed
arm must beat the control by **>~0.5 MAE**.

- Best arm `c64-w128`: 8.117 → 7.685 = **0.432** — **does not clear it.**
- `c64-w128` is also flanked by worse results on both sides (w64 8.117, w256
  8.478), i.e. non-monotonic, which is what a lucky draw looks like rather than
  a trend.
- It additionally *lost* on sex (0.827 vs control 0.865), so it is not even a
  clean within-arm win.

## Conclusion

1. **Keep the locked N-light shape** (`cpg_hidden`/`phi`/`rho` = 64). There is
   no evidence that scaling it up helps once CpGPT is present, and every wider
   arm costs ~15–20% more wall time per epoch.
2. **The G2 CpGPT result stands on its own architecture.** CpGPT's age/sex gain
   was not an artifact of a starved encoder, and is not enlarged by relieving
   it — so the G2 verdict needs no revision.
3. **"Just make it wider" is closed** as a direction for this encoder at this
   scale. Remaining levers are data/representation (full column width,
   non-coding CpGs via cascade), not model size.

## Limits

- Single seed per arm; the design resolves only >~0.5 MAE effects. A *small*
  real capacity benefit (<0.5 MAE) would be invisible here and is not excluded.
- `cpg_hidden` capped at 128 (256 would force batch 64 and ~4× wall time once
  matched), so very wide per-CpG encoders are untested.
- Fold 0, 65k-prefix, N-light only. Says nothing about cascade or full width.
- v1 of this sweep (batch 256, `auto`) was abandoned — see `README.md` for the
  OOM and the batch-size confound that invalidated it.
