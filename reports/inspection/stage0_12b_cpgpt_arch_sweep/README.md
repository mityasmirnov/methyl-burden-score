# CpGPT architecture sweep

## v1 — ABANDONED after 2/6 variants (OOM + batch-size confound)

`ledger_v1_batch256_ABANDONED.json`. Two independent problems:

1. **OOM.** `w128` (cpg_hidden 128) died at batch 256 on the 49 GB GPU0:
   *"Tried to allocate 6.21 GiB … 11.03 GiB is reserved by PyTorch but
   unallocated"*. The chain never set
   `PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True`, which **five** existing
   project scripts set as convention (`run_7h_ats_pooling_s2.sh`,
   `run_7h_next_queue.sh`, `run_7h_dense_stage1_queue.sh`,
   `run_7h_overnight_extension.sh`, `supervise_7g_prime_weekend.py`). Omission
   on our side, and the error message itself named the fix.
2. **Batch-size confound (the more serious one).** v1 used
   `batch_size: auto`. Peak activation memory scales ~ `n_edges (57 430) × batch
   × cpg_hidden`, so wider variants calibrate to *smaller* batches — meaning any
   measured difference would conflate capacity with batch size. This would have
   quietly invalidated the comparison even where it did not OOM.

v1's one completed run (`c64`-equivalent at batch 256, nested age **8.482**) is
retained only as the 7th draw used to correct the run-to-run spread estimate in
`../stage0_12b_cpgpt_multirestart/analysis.md` (Finding 4). It is **not**
comparable to v2 (different batch).

## v2 — current

- `PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True`.
- **`batch_size` fixed at 128 for every variant** (explicit, not `auto`), so
  capacity is the only thing that varies. Recorded per row in the ledger.
- `cpg_hidden` capped at 128: it is the memory-critical per-edge dim, and 256
  would force batch 64, roughly quadrupling wall time once matched across arms.
- Variants: `c64-w{64,128,256}` (downstream capacity only) and
  `c128-w{128,256}` + `c128-w256-d3` (widens the per-CpG encoder — the
  152→cpg_hidden compression that motivated the sweep).

**Resolution limit — read before interpreting.** Run-to-run spread for a *fixed*
config is ~0.53 MAE wide (nested age; see multirestart Finding 4). Only gaps
**>~0.5 MAE** versus the in-sweep control (`c64-w64`) are credible on a single
seed. Smaller orderings are ties, and any apparent winner needs a multi-seed
rerun before promotion.
