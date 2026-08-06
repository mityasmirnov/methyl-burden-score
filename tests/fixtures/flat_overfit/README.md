# Flat baseline overfit fixture

Milestone 5 uses an **in-memory** synthetic fixture from
`mbs.training.dataset.make_synthetic_overfit_bundle` (no on-disk matrix).

```bash
uv run mbs train flat --overfit-fixture --device cpu --max-epochs 200
```
