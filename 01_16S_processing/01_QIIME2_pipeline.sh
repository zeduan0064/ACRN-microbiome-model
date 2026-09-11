#!/bin/bash

###############################################################################
# 16S rRNA V3-V4 sequencing analysis
# QIIME 2 + DADA2 + Greengenes2
#
# Software / database:
#   QIIME 2 v2024.2
#   DADA2
#   Greengenes2 database v2022.10
#
# Notes:
# 1. This script retains the core workflow required for reproducibility.
# 2. All sample identifiers and paths are anonymized/example values.
# 3. No participant identifiers or private server paths are included.
###############################################################################

set -e


###############################################################################
# 1. Activate QIIME 2 environment
###############################################################################

conda activate qiime2-amplicon-2024.2


###############################################################################
# 2. Import paired-end sequencing data
###############################################################################

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.tsv \
  --output-path demux.qza \
  --input-format PairedEndFastqManifestPhred33V2


###############################################################################
# 3. Remove V3-V4 primers
#
# 341F: CCTACGGGNGGCWGCAG
# 805R: GACTACHVGGGTATCTAATCC
###############################################################################

qiime cutadapt trim-paired \
  --i-demultiplexed-sequences demux.qza \
  --p-front-f CCTACGGGNGGCWGCAG \
  --p-front-r GACTACHVGGGTATCTAATCC \
  --p-error-rate 0.1 \
  --p-minimum-length 140 \
  --p-discard-untrimmed \
  --p-cores 30 \
  --o-trimmed-sequences trimmed-seqs.qza


###############################################################################
# 4. Inspect sequencing quality
###############################################################################

qiime demux summarize \
  --i-data trimmed-seqs.qza \
  --o-visualization trimmed-seqs.qzv


###############################################################################
# 5. DADA2 denoising and ASV generation
#
# Truncation parameters should correspond to the sequencing quality profile
# used in the original analysis.
###############################################################################

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


###############################################################################
# 6. Remove low-frequency ASVs occurring in only one sample
###############################################################################

qiime feature-table filter-features \
  --i-table table.qza \
  --p-min-samples 2 \
  --o-filtered-table filtered-table.qza


qiime feature-table filter-seqs \
  --i-data rep-seqs.qza \
  --i-table filtered-table.qza \
  --o-filtered-data filtered-rep-seqs.qza


###############################################################################
# 7. Inspect DADA2 denoising statistics
###############################################################################

qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv


###############################################################################
# 8. Greengenes2 v2022.10 taxonomic annotation
#
# Required Greengenes2 reference files:
#
# refs/2022.10.backbone.full-length.fna.qza
# refs/2022.10.taxonomy.asv.nwk.qza
#
# V3-V4 is treated using the Greengenes2 non-v4-16s workflow.
###############################################################################

qiime greengenes2 non-v4-16s \
  --i-table filtered-table.qza \
  --i-sequences filtered-rep-seqs.qza \
  --i-backbone refs/2022.10.backbone.full-length.fna.qza \
  --o-mapped-table gg2-table.qza \
  --o-representatives gg2-rep-seqs.qza


###############################################################################
# 9. Assign Greengenes2 taxonomy
###############################################################################

qiime greengenes2 taxonomy-from-table \
  --i-reference-taxonomy refs/2022.10.taxonomy.asv.nwk.qza \
  --i-table gg2-table.qza \
  --o-classification taxonomy.qza


###############################################################################
# 10. Visualize taxonomy
###############################################################################

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv


###############################################################################
# 11. Export Greengenes2-mapped feature table
###############################################################################

mkdir -p exported-table

qiime tools export \
  --input-path gg2-table.qza \
  --output-path exported-table

biom convert \
  -i exported-table/feature-table.biom \
  -o exported-table/feature-table.tsv \
  --to-tsv


###############################################################################
# 12. Export taxonomy
###############################################################################

mkdir -p exported-taxonomy

qiime tools export \
  --input-path taxonomy.qza \
  --output-path exported-taxonomy


###############################################################################
# 13. Optional: collapse taxonomy to genus level
#
# Greengenes2 taxonomy levels:
# level 1 = domain
# level 2 = phylum
# level 3 = class
# level 4 = order
# level 5 = family
# level 6 = genus
# level 7 = species
###############################################################################

qiime taxa collapse \
  --i-table gg2-table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table genus-table.qza


###############################################################################
# 14. Export genus-level abundance table
###############################################################################

mkdir -p exported-genus

qiime tools export \
  --input-path genus-table.qza \
  --output-path exported-genus

biom convert \
  -i exported-genus/feature-table.biom \
  -o exported-genus/genus-table.tsv \
  --to-tsv


###############################################################################
# End
###############################################################################

echo "16S rRNA analysis completed."
echo "Taxonomic reference: Greengenes2 v2022.10"
