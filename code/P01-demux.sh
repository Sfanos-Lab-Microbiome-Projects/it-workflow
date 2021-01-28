#!/bin/bash
#SBATCH -c 5
#SBATCH -t 10:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=lpeiffe1@jhu.edu
#SBATCH --job-name=demux atcc-pilot-v2

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
  --input-path  ../data/manifests/chip$i/chip$i-$REGION-manifest.tsv \
  --output-path $outDir/$REGION/chip$i-$REGION-demux.qza \
  --input-format SingleEndFastqManifestPhred33V2

qiime demux summarize \
 --i-data $outDir/$REGION/chip$i-$REGION-demux.qza \
 --o-visualization $outDir/$REGION/chip$i-$REGION-demux.qzv

done
