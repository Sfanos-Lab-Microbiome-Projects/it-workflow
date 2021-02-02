#!/bin/bash
#SBATCH -c 5
#SBATCH -t 24:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=cjone228@jhu.edu
#SBATCH --job-name="P15-clust-tree-cm-goods"
#SBATCH --mem=5G
#SBATCH -p parallel
#SBATCH --array=1-6

# configure
source ~/.bashrc
source activate qiime2-2020.6

outDir=../analysis/P15-clust-tree-cm-goods
mkdir -p $outDir

# array for different regions ${SLURM_ARRAY_TASK_ID} goes from 1 to 6
declare -a arr=("v2" "v3" "v4" "v67" "v8" "v9")
#TASKID=1
REGION=${arr[${SLURM_ARRAY_TASK_ID}-1]}
echo "$REGION"
mkdir -p $outDir/$REGION

# open reference otu clustering
# refer to the P12-dada2 outputs
intableqza=../external/P12-dada2/$REGION/patient-chip1-$REGION-table.qza
inrepsqza=../external/P12-dada2/$REGION/patient-chip1-$REGION-rep-seqs.qza

echo "$intableqza"
echo "$inrepsqza"

qiime vsearch cluster-features-open-reference \
  --i-table $intableqza \
  --i-sequences $inrepsqza \
  --i-reference-sequences ../external/db/sfanos_db_v4.0.fasta.qza \
  --p-perc-identity 0.99 \
  --p-strand plus \
  --p-threads 4 \
  --o-clustered-table $outDir/$REGION/table.qza \
  --o-clustered-sequences $outDir/$REGION/rep-seqs.qza \
  --o-new-reference-sequences $outDir/$REGION/new-ref-seqs.qza

#convert otu table, rep-seqs, and new-ref-seqs from qza to qzv
qiime feature-table summarize \
--i-table ../analysis/P15-clust-tree-cm-goods/$REGION/table.qza \
--o-visualization ../analysis/P15-clust-tree-cm-goods/$REGION/table.qzv \
--m-sample-metadata-file ../data/patient-metadata.tsv \

qiime feature-table tabulate-seqs \
  --i-data ../analysis/P15-clust-tree-cm-goods/$REGION/rep-seqs.qza \
  --o-visualization ../analysis/P15-clust-tree-cm-goods/$REGION/rep-seqs.qzv \
  
qiime feature-table tabulate-seqs \
  --i-data ../analysis/P15-clust-tree-cm-goods/$REGION/new-ref-seqs.qza \
  --o-visualization ../analysis/P15-clust-tree-cm-goods/$REGION/new-ref-seqs.qzv \

mkdir -p ../data/P15-clust-tree-cm-goods/
mkdir -p ../data/P15-clust-tree-cm-goods/$REGION

# make mafft trees with db 4.0
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ../analysis/P15-clust-tree-cm-goods/$REGION/rep-seqs.qza \
  --p-n-threads 4 \
  --o-alignment ../data/P15-clust-tree-cm-goods/$REGION/$REGION-alignment.qza \
  --o-masked-alignment ../data/P15-clust-tree-cm-goods/$REGION/$REGION-masked-alignment.qza \
  --o-tree ../data/P15-clust-tree-cm-goods/$REGION/$REGION-unrooted-tree.qza \
  --o-rooted-tree ../data/P15-clust-tree-cm-goods/$REGION/$REGION-rooted-tree.qza \
  --verbose

#core metrics 
qiime diversity core-metrics-phylogenetic \
--i-table ../analysis/P15-clust-tree-cm-goods/$REGION/table.qza \
--i-phylogeny ../data/P15-clust-tree-cm-goods/$REGION/$REGION-rooted-tree.qza \
--m-metadata-file ../data/patient-metadata.tsv \
--p-sampling-depth 10000 \
--output-dir ../analysis/P15-clust-tree-cm-goods/$REGION/core-metrics-results \
--verbose

# export alpha diversity files
qiime tools export \
  --input-path $outDir/$REGION/core-metrics-results/faith_pd_vector.qza \
  --output-path $outDir/$REGION/exported-faith-pd-vector

qiime tools export \
  --input-path $outDir/$REGION/core-metrics-results/evenness_vector.qza \
  --output-path $outDir/$REGION/exported-evenness-vector
  
qiime tools export \
  --input-path $outDir/$REGION/core-metrics-results/shannon_vector.qza \
  --output-path $outDir/$REGION/exported-shannon-vector
  
qiime tools export \
  --input-path $outDir/$REGION/core-metrics-results/observed_otus_vector.qza \
  --output-path $outDir/$REGION/exported-observed-otus-vector

#taxonomic classification
qiime feature-classifier classify-consensus-vsearch \
  --i-query ../analysis/P15-clust-tree-cm-goods/$REGION/rep-seqs.qza \
  --i-reference-reads ../external/db/sfanos_db_v4.0.fasta.qza \
  --i-reference-taxonomy ../external/db/sfanos_db_v4.0.txt.qza \
  --p-perc-identity 0.99 \
  --o-classification ../analysis/P15-clust-tree-cm-goods/$REGION/tax-class.qza \
  --verbose

# export taxonomic classification qza files
qiime tools export \
  --input-path $outDir/$REGION/tax-class.qza \
  --output-path $outDir/$REGION/tax-class

# calculate Good's Coverage
qiime diversity alpha \
  --i-table $outDir/$REGION/core-metrics-results/rarefied_table.qza \
  --p-metric goods_coverage \
  --o-alpha-diversity $outDir/$REGION/observed_otus_vector.qza
  
# export Good's Coverage qza files
qiime tools export \
  --input-path $outDir/$REGION/observed_otus_vector.qza \
  --output-path $outDir/$REGION/exported-goods-coverage
