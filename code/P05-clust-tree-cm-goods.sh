#!/bin/bash
#SBATCH -c 5
#SBATCH -t 24:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=cjone228@jhu.edu
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
  --i-table ../analysis/P05-clust-tree-cm-goods/$REGION/table.qza \
  --i-taxonomy ../analysis/P05-clust-tree-cm-goods/$REGION/tax-class.qza \
  --p-mode exact \
  --p-exclude "k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales; f__Actinomycetaceae; g__Varibaculum,k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Corynebacteriales; f__Corynebacteriaceae; g__Corynebacterium_1; s__Corynebacterium_amycolatum,k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Micrococcales; f__Microbacteriaceae; g__Microbacterium,k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides; s__Bacteroides_dorei,k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides; s__Bacteroides_faecis,k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides; s__Bacteroides_fragilis,k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides; s__Bacteroides_massiliensis,k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Porphyromonadaceae; g__Parabacteroides; s__Parabacteroides_merdae,k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Bacillaceae; g__Bacillus,k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Family_XI; g__Gemella; s__Gemella_morbillorum,k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Family_XI; g__Anaerococcus; s__Anaerococcus_vaginalis,k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Family_XI; g__Ezakiella; s__Sporobacterium_sp,k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Family_XI; g__Finegoldia; s__Finegoldia_magna,k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Family_XI; g__Peptoniphilus; s__Candidatus_Peptoniphilus_massiliensis,k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Anaerostipes; s__Anaerostipes_hadrus,k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Blautia,k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Lachnoclostridium,k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae,k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Ruminococcaceae; g__Ruminococcus_1; s__Ruminococcus_bicirculans,k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Methylobacteriaceae; g__Methylobacterium,k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Cupriavidus; s__Cupriavidus_metallidurans,k__Bacteria; p__Proteobacteria; c__Deltaproteobacteria; o__Desulfovibrionales; f__Desulfovibrionaceae; g__Bilophila; s__Bilophila_wadsworthia,k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Neisseriales; f__Neisseriaceae; g__Uruburuella; s__Uruburuella_sp,k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacteriales; f__Enterobacteriaceae; g__Cronobacter; s__Cronobacter_sakazakii" \
  --o-filtered-table ../analysis/P05-clust-tree-cm-goods/$REGION/$REGION-filtered-table.qza \

qiime feature-table summarize \
  --i-table ../analysis/P05-clust-tree-cm-goods/$REGION/$REGION-filtered-table.qza \
  --o-visualization ../analysis/P05-clust-tree-cm-goods/$REGION/$REGION-filtered-table.qzv \
  --m-sample-metadata-file ../data/atcc-metadata.tsv \
  
#core metrics 
qiime diversity core-metrics-phylogenetic \
--i-table ../analysis/P05-clust-tree-cm-goods/$REGION/$REGION-filtered-table.qza \
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
