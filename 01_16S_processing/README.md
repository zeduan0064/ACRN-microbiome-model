# 16S rRNA sequencing analysis

本目录包含本研究 16S rRNA 基因 V3-V4 区域测序数据的主要分析代码。

## Software

- QIIME 2 v2024.2
- DADA2
- Greengenes2 database v2022.10

## Workflow

1. 双端 FASTQ 数据导入
2. 引物去除
3. 测序质量评估
4. DADA2 去噪及 ASV 构建
5. 低频 ASV 过滤
6. 使用 Greengenes2 v2022.10 进行物种分类注释
7. 导出 ASV 丰度表
8. 导出 taxonomy 注释结果

## Sequencing region

细菌 16S rRNA 基因 V3-V4 区域使用以下引物扩增：

- 341F: `5'-CCTACGGGNGGCWGCAG-3'`
- 805R: `5'-GACTACHVGGGTATCTAATCC-3'`

测序采用双端测序模式。

## Files

- `01_QIIME2_pipeline.sh`：16S rRNA 数据预处理、DADA2 去噪及 Greengenes2 物种注释流程
- `manifest_example.tsv`：匿名化 manifest 示例文件
- `metadata_example.tsv`：匿名化 metadata 示例文件

## Privacy

为保护研究参与者隐私，本代码库不包含：

- 真实样本编号
- 受试者个人信息
- 原始测序文件
- 临床数据
- 真实服务器路径
- 其他可识别研究参与者身份的信息

示例文件中的样本编号均为虚拟编号，仅用于说明数据格式。

## Notes

DADA2 截断参数应根据实际测序质量结果确定。

物种分类注释使用 Greengenes2 database v2022.10。

如使用本代码分析其他数据集，应根据相应测序质量和研究设计重新确定相关参数。

公开代码仅保留论文复现所需的核心分析流程。
