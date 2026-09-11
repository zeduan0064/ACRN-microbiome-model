# Six machine-learning algorithms

本目录包含最终47菌属面板的六种机器学习算法比较。

## Algorithms

1. Logistic regression
2. Neural network
3. XGBoost
4. Gradient boosting machine
5. Support vector machine
6. AdaBoost

## Core analysis

- repeated 10-fold cross-validation
- 5 repeats
- 1–47 genera incremental AUC comparison
- final 47-genera ROC/AUC overlay
- final 47-genera confusion matrices
- sensitivity, specificity, precision and F1 score
- Logistic regression versus the other five models using DeLong test

## Input

`data/model_input_47genera.csv`

要求：

- 47个菌属列必须按照研究中预设的加入顺序排列
- `Result` 为 0/1 或 No/Yes
- 可包含匿名 `SampleID`

## Output

`results/machine_learning/`

主要结果：

- `Incremental_AUC_6Models_AllResults.csv`
- `Incremental_AUC_6Models.pdf`
- `Incremental_AUC_6Models_CV.pdf`
- `Final_47genera_predictions.csv`
- `Final_47genera_6Models_performance.csv`
- `Final_47genera_ROC_AUC_overlay.pdf`
- `Final_47genera_6Models_ConfusionMatrices.pdf`
- `Final_47genera_confusion_matrix_metrics.csv`
- `Logistic_vs_other_models_DeLong.csv`

## Privacy

公开版本不包含真实样本编号、患者信息、本地路径或服务器路径。