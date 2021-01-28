#!/bin/bash
#SBATCH -c 4
#SBATCH -t 10:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=jwhit117@jhu.edu
#SBATCH --job-name=vsearch-tx-assign
#SBATCH --mem=3G
#SBATCH -p parallel
#SBATCH --account=ksandel1

# configure
source ~/.bashrc
source activate qiime2-2020.2

outDir=../analysis/P00-tx-validation
rm -r $outDir
mkdir -p $outDir

# create data qza
seqfna=../../primer-simulations/analysis/A04-simulate-seqs/extracted.1prcterr.fasta
#seqfna=../../primer-simulations/analysis/A04-simulate-seqs/test.fasta
inrepsqza=../../primer-simulations/analysis/A04-simulate-seqs/extracted.1prcterr.fasta.qza

qiime tools import \
  --input-path $seqfna \
  --output-path $inrepsqza \
  --type 'FeatureData[Sequence]'

echo "$inrepsqza"

qiime feature-classifier classify-consensus-vsearch \
  --i-query $inrepsqza \
  --i-reference-reads ../../db-devel/sfanos_db_v4.0.fasta.qza \
  --i-reference-taxonomy ../../db-devel/sfanos_db_v4.0.txt.qza \
  --p-perc-identity 0.97 \
  --p-threads 4 \
  --o-classification $outDir/classifications.qza \
  --verbose

qiime tools export \
  --input-path $outDir/classifications.qza \
  --output-path $outDir/classifications
