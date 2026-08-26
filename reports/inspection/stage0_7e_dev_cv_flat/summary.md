# Milestone 7E development CV

Split: `hub-ats-7e-3fold-v1` (3 outer folds).

**Winner for Milestone 7:** `N-flat-gene-l1b-nocpgpt` (tissue macro-F1=0.16911202622574176, age RMSE=0.8797459121517698).

Selection: Among neural arms: highest mean tissue macro-F1 on held-out studies, ties broken by lowest age RMSE. Transparent and metadata-only are ceilings.

## Arms

- `T-mean-gene` fold=0 restart=0 L1=None CpGPT=False: tissue_f1=0.17534452472790898 age_rmse=24.056877282470328
- `T-enet` fold=0 restart=0 L1=None CpGPT=False: tissue_f1=0.2907536169361878 age_rmse=26.06256675996183
- `N-flat-gene-l1a` fold=0 restart=0 L1=False CpGPT=True: tissue_f1=0.07707181392565163 age_rmse=1.1121774906721276
- `N-flat-gene-l1b` fold=0 restart=0 L1=True CpGPT=True: tissue_f1=0.18316101247377056 age_rmse=1.0346375026630268
- `N-flat-gene-l1b-nocpgpt` fold=0 restart=0 L1=True CpGPT=False: tissue_f1=0.12282251412064915 age_rmse=0.9381259802783587
- `N-flat-gene-l1a` fold=0 restart=1 L1=False CpGPT=True: tissue_f1=0.09794730174465556 age_rmse=1.1729184478079937
- `N-flat-gene-l1b` fold=0 restart=1 L1=True CpGPT=True: tissue_f1=0.12632855699101128 age_rmse=1.1206180012552653
- `N-flat-gene-l1b-nocpgpt` fold=0 restart=1 L1=True CpGPT=False: tissue_f1=0.17469867172703152 age_rmse=0.9489079574237167
- `T-mean-gene` fold=1 restart=0 L1=None CpGPT=False: tissue_f1=0.2723520591073106 age_rmse=15.915044358377278
- `T-enet` fold=1 restart=0 L1=None CpGPT=False: tissue_f1=0.3171453927853333 age_rmse=22.42842469539218
- `N-flat-gene-l1a` fold=1 restart=0 L1=False CpGPT=True: tissue_f1=0.0467110451971462 age_rmse=0.9157291254505554
- `N-flat-gene-l1b` fold=1 restart=0 L1=True CpGPT=True: tissue_f1=0.1835690628564851 age_rmse=0.772976577856656
- `N-flat-gene-l1b-nocpgpt` fold=1 restart=0 L1=True CpGPT=False: tissue_f1=0.20713397271333164 age_rmse=0.9009387459251351
- `N-flat-gene-l1a` fold=1 restart=1 L1=False CpGPT=True: tissue_f1=0.019363373821840765 age_rmse=1.0678477951087038
- `N-flat-gene-l1b` fold=1 restart=1 L1=True CpGPT=True: tissue_f1=0.1807316465076587 age_rmse=0.8649963629232458
- `N-flat-gene-l1b-nocpgpt` fold=1 restart=1 L1=True CpGPT=False: tissue_f1=0.1842585337595686 age_rmse=0.6706376891699767
- `T-mean-gene` fold=2 restart=0 L1=None CpGPT=False: tissue_f1=0.2591537095251556 age_rmse=18.66307192224545
- `T-enet` fold=2 restart=0 L1=None CpGPT=False: tissue_f1=0.35688114421450484 age_rmse=24.83765131732006
- `N-flat-gene-l1a` fold=2 restart=0 L1=False CpGPT=True: tissue_f1=0.05184529044475629 age_rmse=1.0574620022560934
- `N-flat-gene-l1b` fold=2 restart=0 L1=True CpGPT=True: tissue_f1=0.037144175096340055 age_rmse=0.9442513180139965
- `N-flat-gene-l1b-nocpgpt` fold=2 restart=0 L1=True CpGPT=False: tissue_f1=0.13277546105813634 age_rmse=0.9022371478678689
- `N-flat-gene-l1a` fold=2 restart=1 L1=False CpGPT=True: tissue_f1=0.037208289534362074 age_rmse=1.119882114871952
- `N-flat-gene-l1b` fold=2 restart=1 L1=True CpGPT=True: tissue_f1=0.11882482798396947 age_rmse=0.9467009572329568
- `N-flat-gene-l1b-nocpgpt` fold=2 restart=1 L1=True CpGPT=False: tissue_f1=0.19298300397573342 age_rmse=0.9176279522455629

## Metadata-only controls

- fold=0: `{"age": {"mae": 8.705376626851058, "rmse": 11.757096143441169, "r2": 0.7931726920260997, "pearson_r": 0.8907227836260008, "spearman_r": 0.8327267026424792}, "tissue": {"macro_f1": 0.5524034672970843, "balanced_accuracy": 0.5531914893617021}}`
- fold=1: `{"age": {"mae": 10.39108786601684, "rmse": 13.19882446351047, "r2": 0.6733045048152382, "pearson_r": 0.8206797876366995, "spearman_r": 0.7491236101553056}, "tissue": {"macro_f1": 0.723404255319149, "balanced_accuracy": 0.723404255319149}}`
- fold=2: `{"age": {"mae": 10.197602781950946, "rmse": 13.981778970073492, "r2": 0.6805570651064187, "pearson_r": 0.8250308286783472, "spearman_r": 0.807440114502475}, "tissue": {"macro_f1": 0.7021089466494508, "balanced_accuracy": 0.7021276595744681}}`
