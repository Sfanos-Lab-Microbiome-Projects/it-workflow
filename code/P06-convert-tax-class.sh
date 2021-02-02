#!/bin/bash
#SBATCH -c 6
#SBATCH -t 4:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=lpeiffe1@jhu.edu
#SBATCH --job-name="P06-convert-tax-class"
#SBATCH --mem=6G
#SBATCH -p lrgmem
#SBATCH --array=1-6

# configure
source ~/.bashrc
source activate qiime2-2020.6

# array for different regions ${SLURM_ARRAY_TASK_ID} goes from 1 to 6
declare -a arr=("v2" "v3" "v4" "v67" "v8" "v9")
#TASKID=1
REGION=${arr[${SLURM_ARRAY_TASK_ID}-1]}
echo "$REGION"
outDir=../analysis/P06-convert-tax-class

qiime tools export \
  --input-path ../analysis/P05-clust-tree-cm-goods/$REGION/table.qza \
  --output-path $outDir/$REGION/exported-feature-table

# we need to reformat the column names in $outDir/$REGION/tax-class/taxonomy.tsv
# before adding to $outDir/$REGION/exported-feature-table/feature-table.biom
perl -pe "s/Feature ID\tTaxon\tConsensus/#OTUID\ttaxonomy\tconfidence/g" ../analysis/P05-clust-tree-cm-goods/$REGION/tax-class/taxonomy.tsv > $outDir/$REGION/exported-feature-table/t.txt
mv $outDir/$REGION/exported-feature-table/t.txt $outDir/$REGION/taxonomy.tsv

biom add-metadata \
-i $outDir/$REGION/exported-feature-table/feature-table.biom \
-o $outDir/$REGION/exported-feature-table/feature-table-w-taxonomy.biom \
--observation-metadata-fp $outDir/$REGION/taxonomy.tsv \
--sc-separated taxonomy

biom convert \
-i $outDir/$REGION/exported-feature-table/feature-table-w-taxonomy.biom \
-o $outDir/$REGION/exported-feature-table/feature-table-w-taxonomy.txt \
--header-key taxonomy \
--to-tsv
