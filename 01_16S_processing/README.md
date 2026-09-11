# 16S rRNA sequencing analysis

本目录包含本研究 16S rRNA 基因 V3-V4 区域测序数据的主要分析代码。

## Software

- QIIME2 2024.2
- DADA2
- SILVA 138

## Workflow

1. 双端 FASTQ 数据导入
2. 引物去除
3. 测序质量评估
4. DADA2 降噪
5. ASV 表构建
6. 低频 ASV 过滤
7. SILVA 数据库物种注释
8. ASV 丰度表及 taxonomy 结果导出

## Files

- `01_QIIME2_pipeline.sh`：16S 数据预处理与物种注释代码
- `manifest_example.tsv`：manifest 示例文件
- `metadata_example.tsv`：metadata 示例文件

## Privacy

为保护研究参与者隐私，本代码库不包含真实样本编号、受试者信息、
原始测序文件、真实服务器路径或其他可识别信息。

示例文件中的样本编号均为虚拟编号。

## Notes

DADA2 截断参数应依据实际测序质量结果确定。
如使用本代码分析其他数据集，应根据相应测序质量重新确定相关参数。

公开代码仅保留论文复现所需的核心分析流程。