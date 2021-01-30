#!/bin/bash
#SBATCH -c 12
#SBATCH -t 24:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=cjone228@jhu.edu
#SBATCH --job-name="dada2"
#SBATCH --mem=50G
#SBATCH -p parallel
#SBATCH --array=1-6

# configure
source ~/.bashrc
source activate qiime2-2020.6
LC_ALL=en_US
export LC_ALL

outDir=../external/P02-dada2
mkdir -p $outDir

# array for different regions ${SLURM_ARRAY_TASK_ID} goes from 1 to 6
declare -a arr=("v2" "v3" "v4" "v67" "v8" "v9")
#TASKID=1
REGION=${arr[${SLURM_ARRAY_TASK_ID}-1]}
echo "$REGION"
mkdir -p $outDir/$REGION

for i in {1..3}
do

# denoise pyro
qiime dada2 denoise-pyro \
  --p-trim-left 0 \
  --p-trunc-len 0 \
  --p-max-ee 5 \
  --i-demultiplexed-seqs ../external/P01-demux/$REGION/atcc-chip$i-$REGION-demux.qza \
  --o-table ../external/P02-dada2/$REGION/chip$i-$REGION-table.qza \
  --o-representative-sequences ../external/P02-dada2/$REGION/chip$i-$REGION-rep-seqs.qza \
  --o-denoising-stats ../external/P02-dada2/$REGION/chip$i-$REGION-denoise-stats.qza \
  --p-n-threads 4 \
  --verbose
  
# export to viewable formats --
qiime tools export \
--input-path ../external/P02-dada2/$REGION/chip$i-$REGION-denoise-stats.qza \
--output-path ../external/P02-dada2/$REGION/exported-denoise-stats/chip$i-$REGION-stats.tsv

echo "tabulating DADA2 summary statistics"

# tabulate DADA2 statistics
# put in P02 directory to go along with original un-merged tables and rep-seqs

qiime metadata tabulate \
  --m-input-file ../external/P02-dada2/$REGION/chip$i-$REGION-denoise-stats.qza \
  --o-visualization ../external/P02-dada2/$REGION/chip$i-$REGION-denoise-stats.qzv
  
done

# Merging DADA2 outputs from different sequencing runs
# but for the same region

echo "Performing table merging"
mtableqza=$outDir/$REGION/merged-$REGION-table.qza
mrepsqza=$outDir/$REGION/merged-$REGION-rep-seqs.qza

qiime feature-table merge \
  --i-tables ../external/P02-dada2/$REGION/chip1-$REGION-table.qza \
  --i-tables ../external/P02-dada2/$REGION/chip2-$REGION-table.qza \
  --i-tables ../external/P02-dada2/$REGION/chip3-$REGION-table.qza \
  --o-merged-table $mtableqza

qiime feature-table merge-seqs \
  --i-data ../external/P02-dada2/$REGION/chip1-$REGION-rep-seqs.qza \
  --i-data ../external/P02-dada2/$REGION/chip2-$REGION-rep-seqs.qza \
  --i-data ../external/P02-dada2/$REGION/chip3-$REGION-rep-seqs.qza \
  --o-merged-data $mrepsqza

# export to viewable formats --
qiime tools export \
--input-path $mtableqza \
--output-path $outDir/$REGION/merged

qiime tools export \
--input-path $mrepsqza \
--output-path $outDir/$REGION/merged
