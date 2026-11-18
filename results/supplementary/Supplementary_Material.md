# Supplementary Material

**Network methods find druggable proteins, not disease genes:
a controlled benchmark across 34 diseases**

## Supplementary Table S1. Disease panel

Thirty-four diseases, pre-specified before any analysis. Identifiers were resolved programmatically through the Open Targets search API. Positive labels are proteins with non-zero `clinical` evidence for that disease that are present in the interaction network.

| Disease | Open Targets ID | Matched term | Associated targets | In network | Positive labels | With genetic evidence |
| --- | --- | --- | --- | --- | --- | --- |
| Alzheimer's disease | MONDO_0004975 | Alzheimer disease | 9964 | 10517 | 311 | 1770 |
| Crohn's disease | MONDO_0005011 | Crohn disease | 6314 |  4751 | 137 |  790 |
| Parkinson's disease | MONDO_0005180 | Parkinson disease | 8628 |  7589 | 276 |  504 |
| age-related macular degeneration | MONDO_0005150 | age-related macular degeneration | 3963 |  3594 | 201 |  241 |
| allergic rhinitis | MONDO_0011786 | allergic rhinitis | 3752 |  3461 |  59 |  714 |
| asthma | MONDO_0004979 | asthma | 7491 |  6129 | 249 | 1424 |
| atopic dermatitis | MONDO_0004980 | atopic eczema | 4461 |  4040 | 122 |  747 |
| atrial fibrillation | MONDO_0004981 | atrial fibrillation | 6359 |  5723 | 169 | 1923 |
| bipolar disorder | MONDO_0004985 | bipolar disorder | 3345 |  2990 | 231 |  791 |
| chronic kidney disease | MONDO_0005300 | chronic kidney disease | 6475 |  5761 | 326 |  864 |
| chronic obstructive pulmonary disease | MONDO_0005002 | chronic obstructive pulmonary disease | 7751 |  6896 | 232 |  834 |
| coronary artery disease | MONDO_0005010 | coronary artery disorder | 6546 |  5775 | 274 | 1646 |
| epilepsy | MONDO_0005027 | epilepsy | 7494 |  6774 | 252 |  624 |
| glaucoma | MONDO_0005041 | glaucoma | 5764 |  5093 | 145 | 1330 |
| gout | MONDO_0005393 | gout | 2281 |  2038 |  92 |  918 |
| heart failure | MONDO_0005252 | heart failure | 5479 |  5032 | 277 |  786 |
| hypercholesterolemia | HP_0003124 | Hypercholesterolemia | 2493 |  2282 |  42 |  966 |
| hypertension | MONDO_0005044 | hypertensive disorder | 9863 |  8655 | 353 | 3503 |
| inflammatory bowel disease | MONDO_0005265 | inflammatory bowel disease | 8529 |  7646 | 294 | 1672 |
| major depressive disorder | MONDO_0002009 | major depressive disorder | 6385 |  5856 | 329 | 1037 |
| migraine | MONDO_0005277 | migraine disorder | 2933 |  2709 | 148 |  526 |
| non-alcoholic fatty liver disease | EFO_1001248 | non-alcoholic fatty liver | 5597 |  5162 | 193 |  329 |
| obesity | MONDO_0011122 | obesity disorder | 9855 |  8737 | 297 | 1953 |
| osteoarthritis | MONDO_0005178 | osteoarthritis | 6043 |  5343 | 246 | 1533 |
| osteoporosis | MONDO_0005298 | osteoporosis | 3948 |  3649 |  33 |  209 |
| psoriasis | MONDO_0005083 | psoriasis | 7138 |  5865 | 225 |  797 |
| rheumatoid arthritis | MONDO_0008383 | rheumatoid arthritis | 8249 |  7251 | 344 | 1607 |
| schizophrenia | MONDO_0005090 | schizophrenia | 7739 |  6690 | 266 | 1525 |
| stroke | MONDO_0005098 | stroke disorder | 7871 |  6876 | 374 | 2623 |
| systemic lupus erythematosus | MONDO_0007915 | systemic lupus erythematosus | 6361 |  5468 | 180 |  945 |
| type 1 diabetes mellitus | MONDO_0005147 | type 1 diabetes mellitus | 5884 |  5306 | 290 |  811 |
| type 2 diabetes mellitus | MONDO_0005148 | type 2 diabetes mellitus | 9798 |  8544 | 328 | 3097 |
| ulcerative colitis | MONDO_0005101 | ulcerative colitis | 8900 |  6960 | 232 | 1034 |

## Supplementary Table S2. Benchmark results by method and evaluation scheme

Mean and standard deviation across 34 diseases x 3 repeats x 5 folds.

| Evaluation | Method | AUROC_mean | AUPRC_mean | Precision@50_mean | Precision@48_mean | Precision@47_mean | Precision@34_mean | Precision@33_mean | Precision@27_mean | Precision@26_mean | AUROC_sd | AUPRC_sd | Precision@50_sd | Precision@48_sd | Precision@47_sd | Precision@34_sd | Precision@33_sd | Precision@27_sd | Precision@26_sd |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| all_nodes | RWR | 0.9572 | 0.3593 | 0.3631 | NA | NA | NA | NA | NA | NA | 0.0165 | 0.1123 | 0.1249 | NA | NA | NA | NA | NA | NA |
| all_nodes | spectral + XGBoost | 0.9529 | 0.5647 | 0.4918 | NA | NA | NA | NA | NA | NA | 0.0233 | 0.1482 | 0.1894 | NA | NA | NA | NA | NA | NA |
| all_nodes | GAT | 0.8688 | 0.2592 | 0.2741 | NA | NA | NA | NA | NA | NA | 0.0889 | 0.1421 | 0.1455 | NA | NA | NA | NA | NA | NA |
| all_nodes | XGBoost (no graph) | 0.8509 | 0.2386 | 0.2470 | NA | NA | NA | NA | NA | NA | 0.0382 | 0.1126 | 0.1185 | NA | NA | NA | NA | NA | NA |
| all_nodes | MNC (CytoHubba) | 0.7589 | 0.0420 | 0.0234 | NA | NA | NA | NA | NA | NA | 0.0532 | 0.0334 | 0.0534 | NA | NA | NA | NA | NA | NA |
| all_nodes | Degree | 0.7585 | 0.0420 | 0.0243 | NA | NA | NA | NA | NA | NA | 0.0532 | 0.0333 | 0.0535 | NA | NA | NA | NA | NA | NA |
| all_nodes | Betweenness | 0.7058 | 0.0354 | 0.0454 | NA | NA | NA | NA | NA | NA | 0.0480 | 0.0157 | 0.0331 | NA | NA | NA | NA | NA | NA |
| degree_matched | spectral + XGBoost | 0.9330 | 0.8516 | NA | NA | NA | NA | NA | NA | NA | 0.0312 | 0.0683 | NA | NA | NA | NA | NA | NA | NA |
| degree_matched | RWR | 0.9097 | 0.7863 | NA | NA | NA | NA | NA | NA | NA | 0.0388 | 0.0861 | NA | NA | NA | NA | NA | NA | NA |
| degree_matched | GAT | 0.8137 | 0.6523 | NA | NA | NA | NA | NA | NA | NA | 0.1137 | 0.1656 | NA | NA | NA | NA | NA | NA | NA |
| degree_matched | XGBoost (no graph) | 0.7798 | 0.5807 | NA | NA | NA | NA | NA | NA | NA | 0.0573 | 0.1022 | NA | NA | NA | NA | NA | NA | NA |
| degree_matched | Betweenness | 0.5260 | 0.2991 | NA | NA | NA | NA | NA | NA | NA | 0.0681 | 0.0783 | NA | NA | NA | NA | NA | NA | NA |
| degree_matched | MNC (CytoHubba) | 0.5027 | 0.2723 | NA | NA | NA | NA | NA | NA | NA | 0.0549 | 0.0416 | NA | NA | NA | NA | NA | NA | NA |
| degree_matched | Degree | 0.5009 | 0.2715 | NA | NA | NA | NA | NA | NA | NA | 0.0538 | 0.0424 | NA | NA | NA | NA | NA | NA | NA |

## Supplementary Table S3. Degree effect, every disease and method

Percentage of above-chance performance removed by degree matching, calculated as (AUROC_all - AUROC_matched) / (AUROC_all - 0.5).

| Disease | Method | AUROC all nodes | AUROC degree-matched | % explained by degree |
| --- | --- | --- | --- | --- |
| Alzheimer's disease | Betweenness | 0.7158 | 0.5355 |  83.6 |
| Alzheimer's disease | Degree | 0.7682 | 0.4998 | 100.1 |
| Alzheimer's disease | GAT | 0.8734 | 0.8081 |  17.5 |
| Alzheimer's disease | MNC (CytoHubba) | 0.7715 | 0.5026 |  99.0 |
| Alzheimer's disease | RWR | 0.9517 | 0.8964 |  12.2 |
| Alzheimer's disease | XGBoost (no graph) | 0.8446 | 0.7543 |  26.2 |
| Alzheimer's disease | spectral + XGBoost | 0.9494 | 0.9281 |   4.7 |
| Crohn's disease | Betweenness | 0.7535 | 0.5534 |  78.9 |
| Crohn's disease | Degree | 0.8028 | 0.5020 |  99.3 |
| Crohn's disease | GAT | 0.7969 | 0.6798 |  39.4 |
| Crohn's disease | MNC (CytoHubba) | 0.8018 | 0.5042 |  98.6 |
| Crohn's disease | RWR | 0.9389 | 0.8382 |  22.9 |
| Crohn's disease | XGBoost (no graph) | 0.8420 | 0.7059 |  39.8 |
| Crohn's disease | spectral + XGBoost | 0.9137 | 0.8514 |  15.1 |
| Parkinson's disease | Betweenness | 0.6705 | 0.5388 |  77.2 |
| Parkinson's disease | Degree | 0.6957 | 0.5013 |  99.3 |
| Parkinson's disease | GAT | 0.9064 | 0.8944 |   3.0 |
| Parkinson's disease | MNC (CytoHubba) | 0.6917 | 0.4980 | 101.1 |
| Parkinson's disease | RWR | 0.9407 | 0.9080 |   7.4 |
| Parkinson's disease | XGBoost (no graph) | 0.8307 | 0.7930 |  11.4 |
| Parkinson's disease | spectral + XGBoost | 0.9497 | 0.9359 |   3.1 |
| age-related macular degeneration | Betweenness | 0.7252 | 0.4992 | 100.4 |
| age-related macular degeneration | Degree | 0.8067 | 0.4993 | 100.2 |
| age-related macular degeneration | GAT | 0.8610 | 0.8171 |  12.1 |
| age-related macular degeneration | MNC (CytoHubba) | 0.8065 | 0.5029 |  99.1 |
| age-related macular degeneration | RWR | 0.9555 | 0.9018 |  11.8 |
| age-related macular degeneration | XGBoost (no graph) | 0.8701 | 0.7874 |  22.3 |
| age-related macular degeneration | spectral + XGBoost | 0.9443 | 0.9302 |   3.2 |
| allergic rhinitis | Betweenness | 0.7169 | 0.5858 |  60.4 |
| allergic rhinitis | Degree | 0.6843 | 0.4949 | 102.8 |
| allergic rhinitis | GAT | 0.7684 | 0.7279 |  15.1 |
| allergic rhinitis | MNC (CytoHubba) | 0.6779 | 0.4948 | 102.9 |
| allergic rhinitis | RWR | 0.9698 | 0.9312 |   8.2 |
| allergic rhinitis | XGBoost (no graph) | 0.8114 | 0.6546 |  50.4 |
| allergic rhinitis | spectral + XGBoost | 0.9355 | 0.9067 |   6.6 |
| asthma | Betweenness | 0.6871 | 0.5114 |  93.9 |
| asthma | Degree | 0.7555 | 0.4984 | 100.6 |
| asthma | GAT | 0.8440 | 0.7849 |  17.2 |
| asthma | MNC (CytoHubba) | 0.7555 | 0.5001 |  99.9 |
| asthma | RWR | 0.9493 | 0.8936 |  12.4 |
| asthma | XGBoost (no graph) | 0.8305 | 0.7802 |  15.2 |
| asthma | spectral + XGBoost | 0.9424 | 0.9269 |   3.5 |
| atopic dermatitis | Betweenness | 0.6860 | 0.5364 |  80.4 |
| atopic dermatitis | Degree | 0.6915 | 0.4974 | 101.4 |
| atopic dermatitis | GAT | 0.8007 | 0.7916 |   3.0 |
| atopic dermatitis | MNC (CytoHubba) | 0.6978 | 0.5075 |  96.2 |
| atopic dermatitis | RWR | 0.9536 | 0.9188 |   7.7 |
| atopic dermatitis | XGBoost (no graph) | 0.8103 | 0.6751 |  43.6 |
| atopic dermatitis | spectral + XGBoost | 0.9245 | 0.9146 |   2.3 |
| atrial fibrillation | Betweenness | 0.6878 | 0.5092 |  95.1 |
| atrial fibrillation | Degree | 0.7489 | 0.4984 | 100.6 |
| atrial fibrillation | GAT | 0.9202 | 0.8741 |  11.0 |
| atrial fibrillation | MNC (CytoHubba) | 0.7506 | 0.5001 | 100.0 |
| atrial fibrillation | RWR | 0.9764 | 0.9528 |   5.0 |
| atrial fibrillation | XGBoost (no graph) | 0.8704 | 0.8501 |   5.5 |
| atrial fibrillation | spectral + XGBoost | 0.9755 | 0.9687 |   1.4 |
| bipolar disorder | Betweenness | 0.6744 | 0.4955 | 102.6 |
| bipolar disorder | Degree | 0.7454 | 0.4972 | 101.2 |
| bipolar disorder | GAT | 0.9354 | 0.8731 |  14.3 |
| bipolar disorder | MNC (CytoHubba) | 0.7480 | 0.5010 |  99.6 |
| bipolar disorder | RWR | 0.9706 | 0.9422 |   6.0 |
| bipolar disorder | XGBoost (no graph) | 0.8861 | 0.8284 |  15.0 |
| bipolar disorder | spectral + XGBoost | 0.9742 | 0.9562 |   3.8 |
| chronic kidney disease | Betweenness | 0.7118 | 0.5274 |  87.1 |
| chronic kidney disease | Degree | 0.7690 | 0.4991 | 100.3 |
| chronic kidney disease | GAT | 0.8790 | 0.7593 |  31.6 |
| chronic kidney disease | MNC (CytoHubba) | 0.7695 | 0.5009 |  99.7 |
| chronic kidney disease | RWR | 0.9468 | 0.8855 |  13.7 |
| chronic kidney disease | XGBoost (no graph) | 0.8634 | 0.7861 |  21.3 |
| chronic kidney disease | spectral + XGBoost | 0.9529 | 0.9254 |   6.1 |
| chronic obstructive pulmonary disease | Betweenness | 0.6810 | 0.4862 | 107.6 |
| chronic obstructive pulmonary disease | Degree | 0.7592 | 0.5005 |  99.8 |
| chronic obstructive pulmonary disease | GAT | 0.8860 | 0.8083 |  20.1 |
| chronic obstructive pulmonary disease | MNC (CytoHubba) | 0.7614 | 0.5059 |  97.7 |
| chronic obstructive pulmonary disease | RWR | 0.9585 | 0.9117 |  10.2 |
| chronic obstructive pulmonary disease | XGBoost (no graph) | 0.8626 | 0.8149 |  13.2 |
| chronic obstructive pulmonary disease | spectral + XGBoost | 0.9526 | 0.9331 |   4.3 |
| coronary artery disease | Betweenness | 0.6981 | 0.5276 |  86.1 |
| coronary artery disease | Degree | 0.7459 | 0.5000 | 100.0 |
| coronary artery disease | GAT | 0.8952 | 0.8168 |  19.8 |
| coronary artery disease | MNC (CytoHubba) | 0.7482 | 0.5043 |  98.3 |
| coronary artery disease | RWR | 0.9643 | 0.9366 |   6.0 |
| coronary artery disease | XGBoost (no graph) | 0.8501 | 0.7734 |  21.9 |
| coronary artery disease | spectral + XGBoost | 0.9647 | 0.9536 |   2.4 |
| epilepsy | Betweenness | 0.7081 | 0.5169 |  91.9 |
| epilepsy | Degree | 0.7786 | 0.5594 |  78.7 |
| epilepsy | GAT | 0.8639 | 0.8354 |   7.8 |
| epilepsy | MNC (CytoHubba) | 0.7798 | 0.5614 |  78.1 |
| epilepsy | RWR | 0.9695 | 0.8990 |  15.0 |
| epilepsy | XGBoost (no graph) | 0.8795 | 0.8534 |   6.9 |
| epilepsy | spectral + XGBoost | 0.9698 | 0.9600 |   2.1 |
| glaucoma | Betweenness | 0.6813 | 0.4793 | 111.4 |
| glaucoma | Degree | 0.7470 | 0.5011 |  99.5 |
| glaucoma | GAT | 0.8276 | 0.8297 |  -0.6 |
| glaucoma | MNC (CytoHubba) | 0.7440 | 0.5003 |  99.9 |
| glaucoma | RWR | 0.9697 | 0.9471 |   4.8 |
| glaucoma | XGBoost (no graph) | 0.8266 | 0.7737 |  16.2 |
| glaucoma | spectral + XGBoost | 0.9538 | 0.9427 |   2.4 |
| gout | Betweenness | 0.7159 | 0.4274 | 133.6 |
| gout | Degree | 0.8469 | 0.4990 | 100.3 |
| gout | GAT | 0.9057 | 0.8624 |  10.7 |
| gout | MNC (CytoHubba) | 0.8444 | 0.5013 |  99.6 |
| gout | RWR | 0.9717 | 0.9265 |   9.6 |
| gout | XGBoost (no graph) | 0.8576 | 0.8587 |  -0.3 |
| gout | spectral + XGBoost | 0.9726 | 0.9450 |   5.8 |
| heart failure | Betweenness | 0.6947 | 0.5240 |  87.7 |
| heart failure | Degree | 0.7316 | 0.4997 | 100.1 |
| heart failure | GAT | 0.8654 | 0.8558 |   2.6 |
| heart failure | MNC (CytoHubba) | 0.7304 | 0.4988 | 100.5 |
| heart failure | RWR | 0.9504 | 0.9100 |   9.0 |
| heart failure | XGBoost (no graph) | 0.8289 | 0.7698 |  18.0 |
| heart failure | spectral + XGBoost | 0.9500 | 0.9377 |   2.8 |
| hypercholesterolemia | Betweenness | 0.8094 | 0.6723 |  44.3 |
| hypercholesterolemia | Degree | 0.7858 | 0.4994 | 100.2 |
| hypercholesterolemia | GAT | 0.8950 | 0.8612 |   8.6 |
| hypercholesterolemia | MNC (CytoHubba) | 0.7833 | 0.4934 | 102.3 |
| hypercholesterolemia | RWR | 0.9746 | 0.9363 |   8.1 |
| hypercholesterolemia | XGBoost (no graph) | 0.8301 | 0.7498 |  24.3 |
| hypercholesterolemia | spectral + XGBoost | 0.9586 | 0.9387 |   4.3 |
| hypertension | Betweenness | 0.7113 | 0.5491 |  76.8 |
| hypertension | Degree | 0.7395 | 0.5023 |  99.1 |
| hypertension | GAT | 0.8711 | 0.8073 |  17.2 |
| hypertension | MNC (CytoHubba) | 0.7376 | 0.5023 |  99.0 |
| hypertension | RWR | 0.9435 | 0.8885 |  12.4 |
| hypertension | XGBoost (no graph) | 0.8357 | 0.7642 |  21.3 |
| hypertension | spectral + XGBoost | 0.9446 | 0.9223 |   5.0 |
| inflammatory bowel disease | Betweenness | 0.7147 | 0.5122 |  94.3 |
| inflammatory bowel disease | Degree | 0.7882 | 0.4996 | 100.1 |
| inflammatory bowel disease | GAT | 0.8531 | 0.7146 |  39.2 |
| inflammatory bowel disease | MNC (CytoHubba) | 0.7891 | 0.5010 |  99.7 |
| inflammatory bowel disease | RWR | 0.9402 | 0.8662 |  16.8 |
| inflammatory bowel disease | XGBoost (no graph) | 0.8449 | 0.7507 |  27.3 |
| inflammatory bowel disease | spectral + XGBoost | 0.9378 | 0.9098 |   6.4 |
| major depressive disorder | Betweenness | 0.6585 | 0.5076 |  95.2 |
| major depressive disorder | Degree | 0.7058 | 0.4995 | 100.3 |
| major depressive disorder | GAT | 0.8541 | 0.8612 |  -2.0 |
| major depressive disorder | MNC (CytoHubba) | 0.7052 | 0.5027 |  98.7 |
| major depressive disorder | RWR | 0.9543 | 0.9299 |   5.4 |
| major depressive disorder | XGBoost (no graph) | 0.8461 | 0.7848 |  17.7 |
| major depressive disorder | spectral + XGBoost | 0.9614 | 0.9506 |   2.3 |
| migraine | Betweenness | 0.6485 | 0.5299 |  79.9 |
| migraine | Degree | 0.6739 | 0.4984 | 100.9 |
| migraine | GAT | 0.9426 | 0.9202 |   5.1 |
| migraine | MNC (CytoHubba) | 0.6777 | 0.5024 |  98.7 |
| migraine | RWR | 0.9771 | 0.9682 |   1.9 |
| migraine | XGBoost (no graph) | 0.8602 | 0.7837 |  21.2 |
| migraine | spectral + XGBoost | 0.9685 | 0.9620 |   1.4 |
| non-alcoholic fatty liver disease | Betweenness | 0.6811 | 0.4928 | 104.0 |
| non-alcoholic fatty liver disease | Degree | 0.7510 | 0.4965 | 101.4 |
| non-alcoholic fatty liver disease | GAT | 0.8600 | 0.7710 |  24.7 |
| non-alcoholic fatty liver disease | MNC (CytoHubba) | 0.7534 | 0.4982 | 100.7 |
| non-alcoholic fatty liver disease | RWR | 0.9500 | 0.8930 |  12.7 |
| non-alcoholic fatty liver disease | XGBoost (no graph) | 0.8442 | 0.7833 |  17.7 |
| non-alcoholic fatty liver disease | spectral + XGBoost | 0.9501 | 0.9305 |   4.3 |
| obesity | Betweenness | 0.6780 | 0.5160 |  91.0 |
| obesity | Degree | 0.7236 | 0.4991 | 100.4 |
| obesity | GAT | 0.9149 | 0.8048 |  26.5 |
| obesity | MNC (CytoHubba) | 0.7259 | 0.5043 |  98.1 |
| obesity | RWR | 0.9605 | 0.9309 |   6.4 |
| obesity | XGBoost (no graph) | 0.8615 | 0.8171 |  12.3 |
| obesity | spectral + XGBoost | 0.9594 | 0.9418 |   3.8 |
| osteoarthritis | Betweenness | 0.6918 | 0.4920 | 104.2 |
| osteoarthritis | Degree | 0.7656 | 0.4990 | 100.4 |
| osteoarthritis | GAT | 0.8858 | 0.8146 |  18.5 |
| osteoarthritis | MNC (CytoHubba) | 0.7673 | 0.5031 |  98.8 |
| osteoarthritis | RWR | 0.9618 | 0.9205 |   8.9 |
| osteoarthritis | XGBoost (no graph) | 0.8489 | 0.7765 |  20.8 |
| osteoarthritis | spectral + XGBoost | 0.9594 | 0.9448 |   3.2 |
| osteoporosis | Betweenness | 0.7885 | 0.6622 |  43.8 |
| osteoporosis | Degree | 0.8096 | 0.5001 | 100.0 |
| osteoporosis | GAT | 0.8901 | 0.8287 |  15.7 |
| osteoporosis | MNC (CytoHubba) | 0.8101 | 0.4877 | 104.0 |
| osteoporosis | RWR | 0.9735 | 0.9567 |   3.5 |
| osteoporosis | XGBoost (no graph) | 0.8726 | 0.7949 |  20.8 |
| osteoporosis | spectral + XGBoost | 0.9363 | 0.9398 |  -0.8 |
| psoriasis | Betweenness | 0.7312 | 0.5131 |  94.3 |
| psoriasis | Degree | 0.8329 | 0.4995 | 100.2 |
| psoriasis | GAT | 0.8316 | 0.7912 |  12.2 |
| psoriasis | MNC (CytoHubba) | 0.8335 | 0.5024 |  99.3 |
| psoriasis | RWR | 0.9539 | 0.8632 |  20.0 |
| psoriasis | XGBoost (no graph) | 0.8829 | 0.7975 |  22.3 |
| psoriasis | spectral + XGBoost | 0.9551 | 0.9165 |   8.5 |
| rheumatoid arthritis | Betweenness | 0.7355 | 0.5048 |  98.0 |
| rheumatoid arthritis | Degree | 0.8079 | 0.4979 | 100.7 |
| rheumatoid arthritis | GAT | 0.8237 | 0.8169 |   2.1 |
| rheumatoid arthritis | MNC (CytoHubba) | 0.8095 | 0.5018 |  99.4 |
| rheumatoid arthritis | RWR | 0.9509 | 0.8650 |  19.1 |
| rheumatoid arthritis | XGBoost (no graph) | 0.8612 | 0.7857 |  20.9 |
| rheumatoid arthritis | spectral + XGBoost | 0.9549 | 0.9250 |   6.6 |
| schizophrenia | Betweenness | 0.6789 | 0.5254 |  85.8 |
| schizophrenia | Degree | 0.7164 | 0.4997 | 100.2 |
| schizophrenia | GAT | 0.9166 | 0.8443 |  17.4 |
| schizophrenia | MNC (CytoHubba) | 0.7193 | 0.5045 |  97.9 |
| schizophrenia | RWR | 0.9538 | 0.9216 |   7.1 |
| schizophrenia | XGBoost (no graph) | 0.8502 | 0.8084 |  11.9 |
| schizophrenia | spectral + XGBoost | 0.9552 | 0.9433 |   2.6 |
| stroke | Betweenness | 0.6891 | 0.5561 |  70.3 |
| stroke | Degree | 0.7168 | 0.4996 | 100.2 |
| stroke | GAT | 0.8869 | 0.8707 |   4.2 |
| stroke | MNC (CytoHubba) | 0.7137 | 0.4991 | 100.4 |
| stroke | RWR | 0.9462 | 0.9138 |   7.2 |
| stroke | XGBoost (no graph) | 0.8464 | 0.7862 |  17.4 |
| stroke | spectral + XGBoost | 0.9525 | 0.9388 |   3.0 |
| systemic lupus erythematosus | Betweenness | 0.7338 | 0.4955 | 101.9 |
| systemic lupus erythematosus | Degree | 0.8316 | 0.4964 | 101.1 |
| systemic lupus erythematosus | GAT | 0.7806 | 0.8073 |  -9.5 |
| systemic lupus erythematosus | MNC (CytoHubba) | 0.8319 | 0.5001 | 100.0 |
| systemic lupus erythematosus | RWR | 0.9690 | 0.9142 |  11.7 |
| systemic lupus erythematosus | XGBoost (no graph) | 0.8879 | 0.7838 |  26.9 |
| systemic lupus erythematosus | spectral + XGBoost | 0.9648 | 0.9297 |   7.6 |
| type 1 diabetes mellitus | Betweenness | 0.7149 | 0.5366 |  83.0 |
| type 1 diabetes mellitus | Degree | 0.7689 | 0.4981 | 100.7 |
| type 1 diabetes mellitus | GAT | 0.8808 | 0.7841 |  25.4 |
| type 1 diabetes mellitus | MNC (CytoHubba) | 0.7699 | 0.5015 |  99.4 |
| type 1 diabetes mellitus | RWR | 0.9458 | 0.8839 |  13.9 |
| type 1 diabetes mellitus | XGBoost (no graph) | 0.8573 | 0.7906 |  18.7 |
| type 1 diabetes mellitus | spectral + XGBoost | 0.9515 | 0.9179 |   7.4 |
| type 2 diabetes mellitus | Betweenness | 0.7131 | 0.5450 |  78.9 |
| type 2 diabetes mellitus | Degree | 0.7521 | 0.4995 | 100.2 |
| type 2 diabetes mellitus | GAT | 0.8526 | 0.7524 |  28.4 |
| type 2 diabetes mellitus | MNC (CytoHubba) | 0.7521 | 0.4990 | 100.4 |
| type 2 diabetes mellitus | RWR | 0.9472 | 0.8943 |  11.8 |
| type 2 diabetes mellitus | XGBoost (no graph) | 0.8378 | 0.7558 |  24.3 |
| type 2 diabetes mellitus | spectral + XGBoost | 0.9540 | 0.9318 |   4.9 |
| ulcerative colitis | Betweenness | 0.7058 | 0.4927 | 103.5 |
| ulcerative colitis | Degree | 0.7836 | 0.4992 | 100.3 |
| ulcerative colitis | GAT | 0.9026 | 0.7825 |  29.8 |
| ulcerative colitis | MNC (CytoHubba) | 0.7848 | 0.5031 |  98.9 |
| ulcerative colitis | RWR | 0.9468 | 0.8752 |  16.0 |
| ulcerative colitis | XGBoost (no graph) | 0.8464 | 0.7613 |  24.6 |
| ulcerative colitis | spectral + XGBoost | 0.9554 | 0.9295 |   5.7 |

## Supplementary Table S4. Contribution of disease-specific information

All values are degree-matched AUROC. The disease-agnostic model was trained only to recognise targets of approved drugs, using topology features, with no disease input of any kind.

| Disease | Positive labels | % labels also general drug targets | Own disease model | No disease information | Another disease's model | Gain from disease information |
| --- | --- | --- | --- | --- | --- | --- |
| osteoporosis |  33 | 97.0 | 0.9398 | 0.8411 | 0.8602 | 0.0987 |
| systemic lupus erythematosus | 180 | 81.7 | 0.9297 | 0.8435 | 0.7063 | 0.0862 |
| psoriasis | 225 | 77.8 | 0.9165 | 0.8402 | 0.8738 | 0.0763 |
| hypercholesterolemia |  42 | 83.3 | 0.9387 | 0.8638 | 0.9072 | 0.0749 |
| gout |  92 | 91.3 | 0.9450 | 0.8756 | 0.9171 | 0.0694 |
| ulcerative colitis | 232 | 84.9 | 0.9295 | 0.8623 | 0.9041 | 0.0672 |
| non-alcoholic fatty liver disease | 193 | 82.4 | 0.9305 | 0.8742 | 0.9131 | 0.0563 |
| age-related macular degeneration | 201 | 92.0 | 0.9302 | 0.8789 | 0.8157 | 0.0513 |
| inflammatory bowel disease | 294 | 85.0 | 0.9098 | 0.8592 | 0.8666 | 0.0506 |
| Crohn's disease | 137 | 83.9 | 0.8514 | 0.8015 | 0.7869 | 0.0499 |
| type 2 diabetes mellitus | 328 | 83.5 | 0.9318 | 0.8867 | 0.8851 | 0.0451 |
| heart failure | 277 | 89.9 | 0.9377 | 0.8940 | 0.8930 | 0.0437 |
| Alzheimer's disease | 311 | 88.4 | 0.9281 | 0.8853 | 0.8852 | 0.0428 |
| coronary artery disease | 274 | 91.6 | 0.9536 | 0.9123 | 0.8704 | 0.0413 |
| hypertension | 353 | 87.5 | 0.9223 | 0.8834 | 0.8787 | 0.0389 |
| type 1 diabetes mellitus | 290 | 84.1 | 0.9179 | 0.8818 | 0.9035 | 0.0361 |
| allergic rhinitis |  59 | 93.2 | 0.9067 | 0.8743 | 0.9209 | 0.0324 |
| Parkinson's disease | 276 | 89.9 | 0.9359 | 0.9037 | 0.8798 | 0.0322 |
| asthma | 249 | 85.9 | 0.9269 | 0.8959 | 0.8719 | 0.0310 |
| rheumatoid arthritis | 344 | 81.7 | 0.9250 | 0.8945 | 0.8241 | 0.0305 |
| major depressive disorder | 329 | 93.9 | 0.9506 | 0.9204 | 0.8722 | 0.0302 |
| atrial fibrillation | 169 | 93.5 | 0.9687 | 0.9397 | 0.9175 | 0.0290 |
| chronic obstructive pulmonary disease | 232 | 85.8 | 0.9331 | 0.9049 | 0.9138 | 0.0282 |
| glaucoma | 145 | 91.0 | 0.9427 | 0.9147 | 0.9328 | 0.0280 |
| stroke | 374 | 94.9 | 0.9388 | 0.9122 | 0.8343 | 0.0266 |
| osteoarthritis | 246 | 88.6 | 0.9448 | 0.9218 | 0.9160 | 0.0230 |
| schizophrenia | 266 | 92.9 | 0.9433 | 0.9203 | 0.9095 | 0.0230 |
| bipolar disorder | 231 | 93.9 | 0.9562 | 0.9334 | 0.9358 | 0.0228 |
| obesity | 297 | 93.3 | 0.9418 | 0.9211 | 0.9025 | 0.0207 |
| chronic kidney disease | 326 | 91.7 | 0.9254 | 0.9068 | 0.8212 | 0.0186 |
| atopic dermatitis | 122 | 83.6 | 0.9146 | 0.9003 | 0.8877 | 0.0143 |
| epilepsy | 252 | 97.2 | 0.9600 | 0.9457 | 0.7322 | 0.0143 |
| migraine | 148 | 98.6 | 0.9620 | 0.9572 | 0.9461 | 0.0048 |

## Supplementary Table S5. Similarity of predictions between diseases

Overlap is the Jaccard index. Prediction overlap should be compared with label overlap, not with the rank correlation, which is a different scale.

| Model | Median rank correlation | Top-100 prediction overlap | Label set overlap | Cross-prediction AUROC |
| --- | --- | --- | --- | --- |
| GAT | 0.632 | 0.114 | 0.27 | 0.828 |
| SPECTRAL | 0.488 | 0.290 | 0.27 | 0.872 |

