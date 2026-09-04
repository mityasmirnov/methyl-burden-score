# Seed-panel audit (fold-0)

- panel_hash: `7aefc81aeb60c9dabd13bf0adf817dab`
- graph_content_hash: `None`
- ok_for_seed_mask_gpu: **False**

## Issues

- graph_content_hash is null/missing (provenance blocker)
- age: stability selected 4096 CpGs == prefilter width 4096 — not demonstrably sparse
- age: strength_cap=10876183559.88889 looks unscaled / nonconverged
- sex: stability selected 4096 CpGs == prefilter width 4096 — not demonstrably sparse
- tissue: stability selected 4096 CpGs == prefilter width 4096 — not demonstrably sparse
- sex_autosome control trait missing

## Per-trait selection

| trait | prefilter | stability seeds | sparsity_ok | strength_cap |
|---|---:|---:|:---:|---:|
| age | 4096 | 4096 | False | 10876183559.88889 |
| sex | 4096 | 4096 | False | 49.667053651458204 |
| tissue | 4096 | 4096 | False | 74.83204129007127 |
