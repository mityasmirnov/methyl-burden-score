# Nine-pack trait adequacy (≥1k per arm)

n_samples=34234

| Trait | Verdict | Detail |
|-------|---------|--------|
| age | PASS | n_labeled=10002 |
| tissue | REVIEW | labeled=13457; levels≥1k=1/64: {'whole blood': 2070} |
| sex | PASS | {'Female': 6238, 'Male': 6207} |
| disease | PASS-naive | true=12218 false=22016 (false≠control) |
| cancer | PASS-naive | true=10101 false=24133 (false≠control) |
| blood | FAIL | true=0 false=34234 (false≠control) |
| brain | FAIL | true=0 false=34234 (false≠control) |

Platform: `{'HM450': 34234}`

Blood/brain masks are currently empty in the phenotype table despite pack membership in the virtual route — pack membership ≠ trait label. Disease/cancer mask_true counts look usable but need a proper case/control definition before GPU.
