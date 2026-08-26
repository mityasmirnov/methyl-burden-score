# Milestone 7E development CV

Split: `hub-ats-7e-3fold-v1` (3 outer folds).

**Winner for Milestone 7:** `N-multipath-l1a` (tissue macro-F1=0.32877947753439823, age RMSE=15.05011242468533).

Selection: Among neural arms: highest mean tissue macro-F1 on held-out studies, ties broken by lowest age RMSE. Transparent and metadata-only are ceilings.

## Arms

- `N-gene-direct-l1a` fold=0 restart=0 L1=False CpGPT=True: tissue_f1=0.25199375038519445 age_rmse=19.53545740216931
- `N-gene-direct-l1b` fold=0 restart=0 L1=True CpGPT=True: tissue_f1=0.25199375038519445 age_rmse=19.535517625692666
- `N-rbs-l1a` fold=0 restart=0 L1=False CpGPT=True: tissue_f1=0.15125413589232867 age_rmse=1.0935465016819959
- `N-rbs-l1b` fold=0 restart=0 L1=True CpGPT=True: tissue_f1=0.2033003441596672 age_rmse=1.014995292776245
- `N-tbs-l1a` fold=0 restart=0 L1=False CpGPT=True: tissue_f1=0.004901653453100866 age_rmse=1.0377347408181283
- `N-tbs-l1b` fold=0 restart=0 L1=True CpGPT=True: tissue_f1=0.028016614575171444 age_rmse=1.087041342288947
- `N-multipath-l1a` fold=0 restart=0 L1=False CpGPT=True: tissue_f1=0.27606315554773964 age_rmse=16.7519246213705
- `N-multipath-l1b` fold=0 restart=0 L1=True CpGPT=True: tissue_f1=0.27606315554773964 age_rmse=16.7519246213705
- `N-gene-direct-l1a` fold=0 restart=1 L1=False CpGPT=True: tissue_f1=0.25199375038519445 age_rmse=19.535517625692666
- `N-gene-direct-l1b` fold=0 restart=1 L1=True CpGPT=True: tissue_f1=0.25199375038519445 age_rmse=19.535517625692666
- `N-rbs-l1a` fold=0 restart=1 L1=False CpGPT=True: tissue_f1=0.1444384629611458 age_rmse=1.061461520979198
- `N-rbs-l1b` fold=0 restart=1 L1=True CpGPT=True: tissue_f1=0.20378238618316655 age_rmse=1.0447902071244153
- `N-tbs-l1a` fold=0 restart=1 L1=False CpGPT=True: tissue_f1=0.01657249395606168 age_rmse=1.104244398975693
- `N-tbs-l1b` fold=0 restart=1 L1=True CpGPT=True: tissue_f1=0.06278646916990333 age_rmse=1.107995309834872
- `N-multipath-l1a` fold=0 restart=1 L1=False CpGPT=True: tissue_f1=0.27606315554773964 age_rmse=16.7519246213705
- `N-multipath-l1b` fold=0 restart=1 L1=True CpGPT=True: tissue_f1=0.27606315554773964 age_rmse=16.7519246213705
- `N-gene-direct-l1a` fold=1 restart=0 L1=False CpGPT=True: tissue_f1=0.26972001700189957 age_rmse=15.859516343599807
- `N-gene-direct-l1b` fold=1 restart=0 L1=True CpGPT=True: tissue_f1=0.26972001700189957 age_rmse=15.859516343599807
- `N-rbs-l1a` fold=1 restart=0 L1=False CpGPT=True: tissue_f1=0.04501827013771759 age_rmse=0.9255456253570101
- `N-rbs-l1b` fold=1 restart=0 L1=True CpGPT=True: tissue_f1=0.2439009845782161 age_rmse=0.7605819421797619
- `N-tbs-l1a` fold=1 restart=0 L1=False CpGPT=True: tissue_f1=0.04460181671000782 age_rmse=1.0254862545062453
- `N-tbs-l1b` fold=1 restart=0 L1=True CpGPT=True: tissue_f1=0.03927624932333598 age_rmse=1.0318561629146723
- `N-multipath-l1a` fold=1 restart=0 L1=False CpGPT=True: tissue_f1=0.32744802426285285 age_rmse=14.17959263432762
- `N-multipath-l1b` fold=1 restart=0 L1=True CpGPT=True: tissue_f1=0.32744802426285285 age_rmse=14.17959263432762
- `N-gene-direct-l1a` fold=1 restart=1 L1=False CpGPT=True: tissue_f1=0.26973346840935325 age_rmse=15.859516343599807
- `N-gene-direct-l1b` fold=1 restart=1 L1=True CpGPT=True: tissue_f1=0.26973346840935325 age_rmse=15.859516343599807
- `N-rbs-l1a` fold=1 restart=1 L1=False CpGPT=True: tissue_f1=0.045206153257382885 age_rmse=1.025224957475613
- `N-rbs-l1b` fold=1 restart=1 L1=True CpGPT=True: tissue_f1=0.2184316428361475 age_rmse=0.8196515983023207
- `N-tbs-l1a` fold=1 restart=1 L1=False CpGPT=True: tissue_f1=0.0016923171255135204 age_rmse=0.9530554670948267
- `N-tbs-l1b` fold=1 restart=1 L1=True CpGPT=True: tissue_f1=0.047838471063510496 age_rmse=0.8050738727904849
- `N-multipath-l1a` fold=1 restart=1 L1=False CpGPT=True: tissue_f1=0.32749725761119064 age_rmse=14.17959263432762
- `N-multipath-l1b` fold=1 restart=1 L1=True CpGPT=True: tissue_f1=0.32749725761119064 age_rmse=14.17959263432762
- `N-gene-direct-l1a` fold=2 restart=0 L1=False CpGPT=True: tissue_f1=0.2886760089334398 age_rmse=18.711841340797573
- `N-gene-direct-l1b` fold=2 restart=0 L1=True CpGPT=True: tissue_f1=0.2886760089334398 age_rmse=18.711841340797573
- `N-rbs-l1a` fold=2 restart=0 L1=False CpGPT=True: tissue_f1=0.14990247573240398 age_rmse=0.9614794466716801
- `N-rbs-l1b` fold=2 restart=0 L1=True CpGPT=True: tissue_f1=0.2216763357511087 age_rmse=0.8983902174920394
- `N-tbs-l1a` fold=2 restart=0 L1=False CpGPT=True: tissue_f1=0.02654403134337437 age_rmse=1.0321371807933901
- `N-tbs-l1b` fold=2 restart=0 L1=True CpGPT=True: tissue_f1=0.057450011074225996 age_rmse=0.9370109199869505
- `N-multipath-l1a` fold=2 restart=0 L1=False CpGPT=True: tissue_f1=0.3828026361184334 age_rmse=14.218857124763273
- `N-multipath-l1b` fold=2 restart=0 L1=True CpGPT=True: tissue_f1=0.3828026361184334 age_rmse=14.218857124763273
- `N-gene-direct-l1a` fold=2 restart=1 L1=False CpGPT=True: tissue_f1=0.2886760089334398 age_rmse=18.711841340797573
- `N-gene-direct-l1b` fold=2 restart=1 L1=True CpGPT=True: tissue_f1=0.2886760089334398 age_rmse=18.711841340797573
- `N-rbs-l1a` fold=2 restart=1 L1=False CpGPT=True: tissue_f1=0.0789643267344011 age_rmse=1.0253528448861111
- `N-rbs-l1b` fold=2 restart=1 L1=True CpGPT=True: tissue_f1=0.21745529117630538 age_rmse=0.9650011733619468
- `N-tbs-l1a` fold=2 restart=1 L1=False CpGPT=True: tissue_f1=0.028654058462310634 age_rmse=1.1268291971221789
- `N-tbs-l1b` fold=2 restart=1 L1=True CpGPT=True: tissue_f1=0.05678234152204958 age_rmse=1.010878634828028
- `N-multipath-l1a` fold=2 restart=1 L1=False CpGPT=True: tissue_f1=0.3828026361184334 age_rmse=14.218782911952479
- `N-multipath-l1b` fold=2 restart=1 L1=True CpGPT=True: tissue_f1=0.3828026361184334 age_rmse=14.218782911952479

## Metadata-only controls

- fold=0: `{"age": {"mae": 8.705376626851058, "rmse": 11.757096143441169, "r2": 0.7931726920260997, "pearson_r": 0.8907227836260008, "spearman_r": 0.8327267026424792}, "tissue": {"macro_f1": 0.5524034672970843, "balanced_accuracy": 0.5531914893617021}}`
- fold=1: `{"age": {"mae": 10.39108786601684, "rmse": 13.19882446351047, "r2": 0.6733045048152382, "pearson_r": 0.8206797876366995, "spearman_r": 0.7491236101553056}, "tissue": {"macro_f1": 0.723404255319149, "balanced_accuracy": 0.723404255319149}}`
- fold=2: `{"age": {"mae": 10.197602781950946, "rmse": 13.981778970073492, "r2": 0.6805570651064187, "pearson_r": 0.8250308286783472, "spearman_r": 0.807440114502475}, "tissue": {"macro_f1": 0.7021089466494508, "balanced_accuracy": 0.7021276595744681}}`
