#!/bin/bash
#SBATCH -c 12
#SBATCH -t 24:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=cjone228@jhu.edu
#SBATCH --job-name="dada2 patient samples"
#SBATCH --mem=50G
#SBATCH -p parallel
#SBATCH --array=1-6

# configure
source ~/.bashrc
source activate qiime2-2020.6
LC_ALL=en_US
export LC_ALL

outDir=../external/P12-dada2
mkdir -p $outDir

# array for different regions ${SLURM_ARRAY_TASK_ID} goes from 1 to 6
declare -a arr=("v2" "v3" "v4" "v67" "v8" "v9")
#TASKID=1
REGION=${arr[${SLURM_ARRAY_TASK_ID}-1]}
echo "$REGION"
mkdir -p $outDir/$REGION

# denoise pyro
qiime dada2 denoise-pyro \
  --p-trim-left 0 \
  --p-trunc-len 0 \
  --p-max-ee 5 \
  --i-demultiplexed-seqs ../external/P11-demux/$REGION/patient-chip1-$REGION-demux.qza \
  --o-table ../external/P12-dada2/$REGION/patient-chip1-$REGION-table.qza \
  --o-representative-sequences ../external/P12-dada2/$REGION/patient-chip1-$REGION-rep-seqs.qza \
  --o-denoising-stats ../external/P12-dada2/$REGION/patient-chip1-$REGION-denoise-stats.qza \
  --p-n-threads 4 \
  --verbose
  
# export to viewable formats --
qiime tools export \
--input-path ../external/P12-dada2/$REGION/patient-chip1-$REGION-denoise-stats.qza \
--output-path ../external/P12-dada2/$REGION/exported-denoise-stats/chip1-$REGION-stats.tsv

echo "tabulating DADA2 summary statistics"

# tabulate DADA2 statistics
# put in P12 directory to go along with tables and rep-seqs

qiime metadata tabulate \
  --m-input-file ../external/P12-dada2/$REGION/patient-chip1-$REGION-denoise-stats.qza \
  --o-visualization ../external/P12-dada2/$REGION/patient-chip1-$REGION-denoise-stats.qzv

# export feature tables and rep seqs to viewable formats --
qiime tools export \
--input-path ../external/P12-dada2/$REGION/patient-chip1-$REGION-table.qza \
--output-path ../external/P12-dada2/$REGION/exported-patient-table

qiime tools export \
--input-path ../external/P12-dada2/$REGION/patient-chip1-$REGION-rep-seqs.qza \
--output-path ..external/P12-dada2/$REGION/exported-patient-rep-seqs
