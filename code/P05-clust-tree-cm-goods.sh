#!/bin/bash
#SBATCH -c 5
#SBATCH -t 24:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=lpeiffe1@jhu.edu
#SBATCH --job-name="P05-clust-tree-cm-goods"
#SBATCH --mem=5G
#SBATCH -p parallel
#SBATCH --array=1-6

# configure
source ~/.bashrc
source activate qiime2-2020.6

outDir=../analysis/P05-clust-tree-cm-goods
mkdir -p $outDir

# array for different regions ${SLURM_ARRAY_TASK_ID} goes from 1 to 6
declare -a arr=("v2" "v3" "v4" "v67" "v8" "v9")
#TASKID=1
REGION=${arr[${SLURM_ARRAY_TASK_ID}-1]}
echo "$REGION"
mkdir -p $outDir/$REGION

# open reference otu clustering
# refer to the P02-dada2 outputs
intableqza=../external/P02-dada2/$REGION/merged-$REGION-table.qza
inrepsqza=../external/P02-dada2/$REGION/merged-$REGION-rep-seqs.qza

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
--i-table ../analysis/P05-clust-tree-cm-goods/$REGION/table.qza \
--o-visualization ../analysis/P05-clust-tree-cm-goods/$REGION/table.qzv \
--m-sample-metadata-file ../data/atcc-metadata.tsv \

qiime feature-table tabulate-seqs \
  --i-data ../analysis/P05-clust-tree-cm-goods/$REGION/rep-seqs.qza \
  --o-visualization ../analysis/P05-clust-tree-cm-goods/$REGION/rep-seqs.qzv \
  
qiime feature-table tabulate-seqs \
  --i-data ../analysis/P05-clust-tree-cm-goods/$REGION/new-ref-seqs.qza \
  --o-visualization ../analysis/P05-clust-tree-cm-goods/$REGION/new-ref-seqs.qzv \

mkdir -p ../data/P05-clust-tree-cm-goods/
mkdir -p ../data/P05-clust-tree-cm-goods/$REGION

# make mafft trees with db 4.0
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ../analysis/P05-clust-tree-cm-goods/$REGION/rep-seqs.qza \
  --p-n-threads 4 \
  --o-alignment ../data/P05-clust-tree-cm-goods/$REGION/$REGION-alignment.qza \
  --o-masked-alignment ../data/P05-clust-tree-cm-goods/$REGION/$REGION-masked-alignment.qza \
  --o-tree ../data/P05-clust-tree-cm-goods/$REGION/$REGION-unrooted-tree.qza \
  --o-rooted-tree ../data/P05-clust-tree-cm-goods/$REGION/$REGION-rooted-tree.qza \
  --verbose

#taxonomic classification
qiime feature-classifier classify-consensus-vsearch \
  --i-query ../analysis/P05-clust-tree-cm-goods/$REGION/rep-seqs.qza \
  --i-reference-reads ../external/db/sfanos_db_v4.0.fasta.qza \
  --i-reference-taxonomy ../external/db/sfanos_db_v4.0.txt.qza \
  --p-perc-identity 0.99 \
  --o-classification ../analysis/P05-clust-tree-cm-goods/$REGION/tax-class.qza \
  --verbose

# export taxonomic classification qza files
qiime tools export \
  --input-path $outDir/$REGION/tax-class.qza \
  --output-path $outDir/$REGION/tax-class

#filter contaminants out of feature table
qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-mode exact \
  --p-exclude "k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales; f__Actinomycetaceae; g__Varibaculum; s__unassigned","s__unassigned, k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Corynebacteriales; f__Corynebacteriaceae; g__Corynebacterium_1; s__Corynebacterium_amycolatum","__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Micrococcales; f__Microbacteriaceae; g__Microbacterium; s__unassigned" \
  --o-filtered-table table-no-mitochondria-exact.qza
  
#core metrics 
qiime diversity core-metrics-phylogenetic \
--i-table ../analysis/P05-clust-tree-cm-goods/$REGION/table.qza \
--i-phylogeny ../data/P05-clust-tree-cm-goods/$REGION/$REGION-rooted-tree.qza \
--m-metadata-file ../data/atcc-metadata.tsv \
--p-sampling-depth 10000 \
--output-dir ../analysis/P05-clust-tree-cm-goods/$REGION/core-metrics-results \
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
  --input-path $outDir/$REGION/core-metrics-results/observed_features_vector.qza \
  --output-path $outDir/$REGION/exported-observed-otus-vector

# export beta diversity files
qiime tools export \
 --input-path $outDir/$REGION/core-metrics-results/bray_curtis_distance_matrix.qza \
 --output-path $outDir/$REGION/core-metrics-results/bray-curtis-dm

qiime tools export \
 --input-path $outDir/$REGION/core-metrics-results/jaccard_distance_matrix.qza \
 --output-path $outDir/$REGION/core-metrics-results/jaccard-dm 

qiime tools export \
 --input-path $outDir/$REGION/core-metrics-results/unweighted_unifrac_distance_matrix.qza \
 --output-path $outDir/$REGION/core-metrics-results/unweighted-unifrac-dm

qiime tools export \
 --input-path $outDir/$REGION/core-metrics-results/weighted_unifrac_distance_matrix.qza \
 --output-path $outDir/$REGION/core-metrics-results/weighted-unifrac-dm

# calculate Good's Coverage
qiime diversity alpha \
  --i-table $outDir/$REGION/core-metrics-results/rarefied_table.qza \
  --p-metric goods_coverage \
  --o-alpha-diversity $outDir/$REGION/observed_otus_vector.qza
  
# export Good's Coverage qza files
qiime tools export \
  --input-path $outDir/$REGION/observed_otus_vector.qza \
  --output-path $outDir/$REGION/goods-coverage
