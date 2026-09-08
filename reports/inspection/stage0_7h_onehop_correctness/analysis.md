# One-hop correctness smokes (Milestone 10)

Generated: `2026-09-08T12:09:01.023864+00:00`

Cheap correctness only. One-hop+seed-mask and one-hop multi-seed were previously untested; do not scale before these land. Uses N-light default rho_hidden=64.

## Seed-mask (fold 0, G0 vs G1 on one-hop)

- **n_seed_genes_age**: `256`
- **n_seed_genes_tissue**: `41`
- **n_seed_genes_sex**: `50`
- **rho_hidden_dim**: `64`
- **G0_onehop**: `{'tissue_f1': 0.3350715881349236, 'age_mae': 17.30601512710046, 'sex_auroc': 0.7503350254464041}`
- **G1_onehop**: `{'tissue_f1': 0.14590206771158448, 'age_mae': 21.712911626753257, 'sex_auroc': 0.5596277351244495}`

### Interpretation

G0 = dense one-hop (no seed mask); G1 = trait-associated seed genes only. On cascade (9c), G0 beat G1–G3 for age-primary training. Here we ask whether one-hop’s single-stage pooling interacts differently with sparse masking. Stable train + sensible metrics = correctness pass; adopt masking only if G1 clearly helps on the primary tasks.

- G0 tissue/age/sex: `{'tissue_f1': 0.3350715881349236, 'age_mae': 17.30601512710046, 'sex_auroc': 0.7503350254464041}`
- G1 tissue/age/sex: `{'tissue_f1': 0.14590206771158448, 'age_mae': 21.712911626753257, 'sex_auroc': 0.5596277351244495}`

## Multi-seed restarts (fold 0)

- **rho_hidden_dim**: `64`
- **onehop_mean_s42**: `{'seed': 42, 'tissue_f1': 0.2988096584577759, 'age_mae': 17.661871715225228, 'sex_auroc': 0.7997804782565088}`
- **onehop_mean_s43**: `{'seed': 43, 'tissue_f1': 0.3254996663074338, 'age_mae': 16.48994905198777, 'sex_auroc': 0.7668924974359471}`
- **onehop_mean_s44**: `{'seed': 44, 'tissue_f1': 0.3114642620650313, 'age_mae': 15.693638692739013, 'sex_auroc': 0.8011877145473194}`
- **_diversity**: `{'tissue_f1_span': 0.026690007849657915, 'n_seeds': 3, 'ok_distinct': True}`

### Interpretation

After the seed-offset fix, restarts must produce **distinct** sensible runs (tissue F1 span `0.026690007849657915`, ok_distinct=`True`). A zero span means the seed path is still broken; do not trust Milestone 12 restart ensembling until this is true.

