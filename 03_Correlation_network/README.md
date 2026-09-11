# Genus-pathway correlation and co-occurrence network

本目录将菌属-通路相关性分析和菌属共现网络分析合并为一个 R 脚本。

## 文件

- `03_correlation_network_analysis.R`

## 输入

脚本默认从项目根目录的 `data/` 文件夹读取：

- `genus_47_relative_abundance.csv`
- `pathway_abundance.csv`
- `metadata.tsv`

### genus_47_relative_abundance.csv

- 行：匿名样本
- 第一列：`SampleID` 或 `sample-id`
- 其余列：目标菌属相对丰度

### pathway_abundance.csv

- 行：Pathway
- 第一列：Pathway名称
- 其余列：匿名样本

### metadata.tsv

至少包含：

- `SampleID` 或 `sample-id`
- `Group`：HC / ACRN

## 分析1：菌属与通路相关性

1. Pathway丰度进行 CLR 转换
2. HC和ACRN分别分析
3. Spearman correlation
4. FDR校正
5. 输出相关系数、FDR和组合热图

## 分析2：菌属共现网络

1. 使用菌属相对丰度
2. prevalence >= 5%
3. Spearman correlation
4. BH-FDR
5. 1000次bootstrap
6. Bootstrap 95% CI不跨0
7. FDR < 0.05
8. |median bootstrap rho| >= 0.20
9. 输出网络边、节点属性、网络统计及网络图

## Privacy

公开版本已删除真实样本编号、真实路径、患者信息以及与结果复现无关的手工绘图参数。
分组信息通过匿名化 metadata 文件提供，不再通过真实样本编号规则推断。