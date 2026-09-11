# Iterative LASSO feature selection

本目录仅包含研究中的迭代 LASSO 特征筛选代码。

## 核心流程

1. 菌属丰度进行 CLR 转换
2. 重复运行 50 次 LASSO
3. 每次随机抽取 70% 样本
4. 每次进行 10-fold cross-validation
5. 使用 `lambda.min` 提取非零系数特征
6. 统计每个菌属在 50 次迭代中的入选次数
7. 按入选频率排序
8. 保留排名前 20% 的候选特征

## 输入

`data/genus_abundance.tsv`

数据格式：

- 第一列：匿名样本编号
- 中间列：菌属丰度
- 最后一列：二分类标签（0/1）

## 输出

结果保存到：

`results/feature_selection/`

主要文件：

- `LASSO_50x_summary.csv`
- `LASSO_top20_percent.csv`

## Privacy

公开版本不包含真实样本编号、患者信息、真实本地路径或服务器路径。