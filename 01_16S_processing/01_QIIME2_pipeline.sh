#!/bin/bash

###############################################################################
# 16S rRNA V3-V4 sequencing analysis
# QIIME2 + DADA2
#
# 说明：
# 1. 本脚本仅保留论文复现所需的核心步骤。
# 2. 所有样本编号和路径均已匿名化/示例化。
# 3. 不应公开真实受试者信息、真实样本编号或服务器绝对路径。
###############################################################################

# 1. 激活 QIIME2 环境
conda activate qiime2-amplicon-2024.2

# 2. 导入双端测序数据
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.tsv \
  --output-path demux.qza \
  --input-format PairedEndFastqManifestPhred33V2

# 3. 去除 V3-V4 引物
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences demux.qza \
  --p-front-f CCTACGGGNGGCWGCAG \
  --p-front-r GACTACHVGGGTATCTAATCC \
  --p-error-rate 0.1 \
  --p-minimum-length 140 \
  --p-discard-untrimmed \
  --p-cores 30 \
  --o-trimmed-sequences trimmed-seqs.qza

# 4. 查看测序质量
qiime demux summarize \
  --i-data trimmed-seqs.qza \
  --o-visualization trimmed-seqs.qzv

# 5. DADA2 去噪并生成 ASV
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs trimmed-seqs.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 240 \
  --p-trunc-len-r 200 \
  --p-n-threads 30 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza

# 6. 去除仅出现在单一样本中的低频 ASV
qiime feature-table filter-features \
  --i-table table.qza \
  --p-min-samples 2 \
  --o-filtered-table filtered-table.qza

qiime feature-table filter-seqs \
  --i-data rep-seqs.qza \
  --i-table filtered-table.qza \
  --o-filtered-data filtered-rep-seqs.qza

# 7. 查看 DADA2 质控结果
qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv

# 8. 导出 ASV 丰度表
mkdir -p exported-table

qiime tools export \
  --input-path filtered-table.qza \
  --output-path exported-table

biom convert \
  -i exported-table/feature-table.biom \
  -o exported-table/feature-table.tsv \
  --to-tsv

# 9. 提取 SILVA V3-V4 参考序列
qiime feature-classifier extract-reads \
  --i-sequences silva-138-99-seqs.qza \
  --p-f-primer CCTACGGGNGGCWGCAG \
  --p-r-primer GACTACHVGGGTATCTAATCC \
  --p-min-length 200 \
  --p-max-length 600 \
  --o-reads silva-v3v4-seqs.qza

# 10. 训练 Naive Bayes 分类器
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-v3v4-seqs.qza \
  --i-reference-taxonomy silva-138-99-tax.qza \
  --o-classifier silva-v3v4-classifier.qza

# 11. 物种注释
qiime feature-classifier classify-sklearn \
  --i-classifier silva-v3v4-classifier.qza \
  --i-reads filtered-rep-seqs.qza \
  --o-classification taxonomy.qza

# 12. 导出物种注释结果
mkdir -p exported-taxonomy

qiime tools export \
  --input-path taxonomy.qza \
  --output-path exported-taxonomy