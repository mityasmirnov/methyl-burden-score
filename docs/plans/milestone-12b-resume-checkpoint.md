# §3.6 Resume-from-checkpoint (flat training loop)

Status: **DONE** 2026-10-07 (flat loop). Cascade resume is out of scope.
Tests: `tests/unit/test_resume_checkpoint.py`.

## Scope and acceptance

**Done when:** flat `train_flat_baseline` can continue a run from `last.pt` after
kill/preemption, and CPU overfit tests prove:

1. Train N straight ≡ train K then resume to N (final weights + `history`).
2. Truncated / corrupt `last.pt` → explicit error (no silent load).
3. `config_hash` mismatch → explicit error.
4. Checkpoint without `resume_state` (pre-§3.6) → explicit error.
5. Failed mid-write leaves the previous `last.pt` intact (atomic replace).

## Locked decisions

| Choice | Decision | Why |
|--------|----------|-----|
| Atomic write | `torch.save` → `*.tmp` then `os.replace` | Crash mid-write must not truncate `last.pt` |
| When to write `last.pt` | After best/patience bookkeeping | `resume_state` matches `best.pt` |
| Resume config | `training.resume: auto \| <path>` | `auto` uses `checkpoints/<run_id>/last.pt` if present |
| Fail closed | Hash / missing `resume_state` / `epoch >= max_epochs` / unreadable file | Repo rule: raise, never guess |
| RNG | Save torch (+ CUDA) and numpy states | Dropout fidelity; no GPU bit-identity promise |
| Calibrate on resume | Skip VRAM probe; restore `batch_size` | Same micro-batch as interrupted run |
| Cascade | Deferred | `cascade_loop.py` has warm-start only |

## `resume_state` schema (inside checkpoint payload)

```text
best_val, best_rank, best_epoch, stale,
history, val_history,
n_samples_seen, n_optimizer_steps,
batch_size, batch_token_budget,
torch_rng_state, numpy_rng_state, cuda_rng_state_all?
```

## Non-goals

- Bit-identical GPU resume (no determinism flags in repo).
- Cascade / hier loop resume.
- Changing `dev_cv.py` skip-finished-run behaviour.

## Open questions

None blocking.
