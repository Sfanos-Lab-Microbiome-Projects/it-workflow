#!/bin/bash
#SBATCH -c 12
#SBATCH -t 10:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=cjone228@jhu.edu
#SBATCH --job-name="demux all atcc"
#SBATCH --mem=20G
#SBATCH -p parallel
#SBATCH --array=1-6

source ~/.bashrc
#source activate qiime2-2020.6

locale
LC_ALL=en_US.UTF-8
export LC_ALL
LANG="en_US.UTF-8"
export LANG
locale

conda activate qiime2-2020.6

outDir=../external/P01-demux
mkdir -p $outDir

# array for different regions ${SLURM_ARRAY_TASK_ID} goes from 1 to 6
declare -a arr=("v2" "v3" "v4" "v67" "v8" "v9")
#TASKID=1
REGION=${arr[${SLURM_ARRAY_TASK_ID}-1]}
echo "$REGION"
mkdir -p $outDir/$REGION

# for loop corresponds to chips
for i in {1..3}
do

qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path  ../data/manifests/atcc-chip$i-$REGION.tsv \
  --output-path $outDir/$REGION/atcc-chip$i-$REGION-demux.qza \
  --input-format SingleEndFastqManifestPhred33V2

qiime demux summarize \
 --i-data $outDir/$REGION/atcc-chip$i-$REGION-demux.qza \
 --o-visualization $outDir/$REGION/atcc-chip$i-$REGION-demux.qzv

done
