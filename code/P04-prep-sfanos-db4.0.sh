#!/bin/bash
#SBATCH -c 1
#SBATCH -t 10:00:00
#SBATCH -p lrgmem
#SBATCH --mail-type=all
#SBATCH --mail-user=cjone228@jhu.edu
#SBATCH --job-name=import-db
#SBATCH --mem=4G
#SBATCH --account=ksandel1

# configure
source ~/.bashrc
source activate qiime2-2020.6
LC_ALL=en_US
export LC_ALL

# new and improved db of reference sequences
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path ../../db-devel/sfanos_db_v4.0.fasta \
  --output-path ../../db-devel/sfanos_db_v4.0.fasta.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-path ../../db-devel/sfanos_db_v4.0.txt \
  --input-format HeaderlessTSVTaxonomyFormat \
  --output-path ../../db-devel/sfanos_db_v4.0.txt.qza
