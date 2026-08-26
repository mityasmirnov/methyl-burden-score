# Milestone 7E development CV

Split: `hub-ats-7e-3fold-v1` (3 outer folds).

**Winner for Milestone 7:** `N-hier-gene-l1b-nocpgpt` (tissue macro-F1=0.22453646614271386, age RMSE=0.9686791202187127).

Selection: Among neural arms: highest mean tissue macro-F1 on held-out studies, ties broken by lowest age RMSE. Transparent and metadata-only are ceilings.

## Arms

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

## Metadata-only controls

- fold=0: `{"age": {"mae": 8.705376626851058, "rmse": 11.757096143441169, "r2": 0.7931726920260997, "pearson_r": 0.8907227836260008, "spearman_r": 0.8327267026424792}, "tissue": {"macro_f1": 0.5524034672970843, "balanced_accuracy": 0.5531914893617021}}`
- fold=1: `{"age": {"mae": 10.39108786601684, "rmse": 13.19882446351047, "r2": 0.6733045048152382, "pearson_r": 0.8206797876366995, "spearman_r": 0.7491236101553056}, "tissue": {"macro_f1": 0.723404255319149, "balanced_accuracy": 0.723404255319149}}`
- fold=2: `{"age": {"mae": 10.197602781950946, "rmse": 13.981778970073492, "r2": 0.6805570651064187, "pearson_r": 0.8250308286783472, "spearman_r": 0.807440114502475}, "tissue": {"macro_f1": 0.7021089466494508, "balanced_accuracy": 0.7021276595744681}}`
