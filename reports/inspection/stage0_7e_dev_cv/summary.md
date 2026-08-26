# Milestone 7E development CV

Split: `hub-ats-7e-3fold-v1` (3 outer folds).
CV budget: max_loci=8192, max_epochs=2.

**Winner for Milestone 7:** `N-multipath-l1a` (tissue macro-F1=0.32877947753439823, age RMSE=15.05011242468533).

Selection: Among neural architecture arms (flat/hier/gene-direct/multipath): highest mean tissue macro-F1 on held-out studies, ties broken by lowest age RMSE. Transparent and metadata-only are ceilings.

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
- `N-hier-gene-l1a` fold=0 restart=0 L1=False CpGPT=True: tissue_f1=0.09229622774444389 age_rmse=1.2036701678140689
- `N-hier-gene-l1b` fold=0 restart=0 L1=True CpGPT=True: tissue_f1=0.1761011440538889 age_rmse=1.0216648125457513
- `N-hier-gene-l1b-nocpgpt` fold=0 restart=0 L1=True CpGPT=False: tissue_f1=0.26417559946882185 age_rmse=0.9927936013642761
- `N-hier-gene-l1a` fold=0 restart=1 L1=False CpGPT=True: tissue_f1=0.06437976689387855 age_rmse=1.2668021445991535
- `N-hier-gene-l1b` fold=0 restart=1 L1=True CpGPT=True: tissue_f1=0.1885145743323245 age_rmse=1.1164484435847493
- `N-hier-gene-l1b-nocpgpt` fold=0 restart=1 L1=True CpGPT=False: tissue_f1=0.22631789654835502 age_rmse=0.998428642362978
- `N-hier-gene-l1a` fold=1 restart=0 L1=False CpGPT=True: tissue_f1=0.13643350408066304 age_rmse=0.9216069561278402
- `N-hier-gene-l1b` fold=1 restart=0 L1=True CpGPT=True: tissue_f1=0.19158966289448712 age_rmse=0.7742890128377308
- `N-hier-gene-l1b-nocpgpt` fold=1 restart=0 L1=True CpGPT=False: tissue_f1=0.2742521992281007 age_rmse=0.8550035473109002
- `N-hier-gene-l1a` fold=1 restart=1 L1=False CpGPT=True: tissue_f1=0.0581433697430244 age_rmse=0.9412569588896379
- `N-hier-gene-l1b` fold=1 restart=1 L1=True CpGPT=True: tissue_f1=0.0707885176779312 age_rmse=0.9922843357581961
- `N-hier-gene-l1b-nocpgpt` fold=1 restart=1 L1=True CpGPT=False: tissue_f1=0.26862436965929914 age_rmse=0.7708062352654977
- `N-hier-gene-l1a` fold=2 restart=0 L1=False CpGPT=True: tissue_f1=0.061874617839102634 age_rmse=1.2213920353661565
- `N-hier-gene-l1b` fold=2 restart=0 L1=True CpGPT=True: tissue_f1=0.19699705670148182 age_rmse=1.003543694318423
- `N-hier-gene-l1b-nocpgpt` fold=2 restart=0 L1=True CpGPT=False: tissue_f1=0.13925493719891308 age_rmse=1.1308505229863004
- `N-hier-gene-l1a` fold=2 restart=1 L1=False CpGPT=True: tissue_f1=0.058951026546788174 age_rmse=1.351620589279349
- `N-hier-gene-l1b` fold=2 restart=1 L1=True CpGPT=True: tissue_f1=0.1388197401969827 age_rmse=1.0695595373573723
- `N-hier-gene-l1b-nocpgpt` fold=2 restart=1 L1=True CpGPT=False: tissue_f1=0.17459379475279327 age_rmse=1.064192172022325
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

- fold=0: `{"age": {"mae": 8.705376626851058, "pearson_r": 0.8907227836260008, "r2": 0.7931726920260997, "rmse": 11.757096143441169, "spearman_r": 0.8327267026424792}, "tissue": {"balanced_accuracy": 0.5531914893617021, "macro_f1": 0.5524034672970843}}`
- fold=1: `{"age": {"mae": 10.39108786601684, "pearson_r": 0.8206797876366995, "r2": 0.6733045048152382, "rmse": 13.19882446351047, "spearman_r": 0.7491236101553056}, "tissue": {"balanced_accuracy": 0.723404255319149, "macro_f1": 0.723404255319149}}`
- fold=2: `{"age": {"mae": 10.197602781950946, "pearson_r": 0.8250308286783472, "r2": 0.6805570651064187, "rmse": 13.981778970073492, "spearman_r": 0.807440114502475}, "tissue": {"balanced_accuracy": 0.7021276595744681, "macro_f1": 0.7021089466494508}}`
