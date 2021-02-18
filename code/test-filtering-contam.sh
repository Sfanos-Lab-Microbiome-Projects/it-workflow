#!/bin/bash
#SBATCH -c 5
#SBATCH -t 24:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=cjone228@jhu.edu
#SBATCH --job-name="test filtering contam"
#SBATCH --mem=5G
#SBATCH -p parallel
#SBATCH --array=1-6

# configure
source ~/.bashrc
source activate qiime2-2020.6

# array for different regions ${SLURM_ARRAY_TASK_ID} goes from 1 to 6
declare -a arr=("v2" "v3" "v4" "v67" "v8" "v9")
#TASKID=1
REGION=${arr[${SLURM_ARRAY_TASK_ID}-1]}
echo "$REGION"
mkdir -p $outDir/$REGION

qiime taxa filter-table \
  --i-table ../analysis/P05-clust-tree-cm-goods/$REGION/table.qza \
  --i-taxonomy ../analysis/P05-clust-tree-cm-goods/$REGION/tax-class.qza \
  --p-mode exact \
  --p-exclude "k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales; f__Actinomycetaceae; g__Varibaculum,k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Corynebacteriales; f__Corynebacteriaceae; g__Corynebacterium_1; s__Corynebacterium_amycolatum,k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Micrococcales; f__Microbacteriaceae; g__Microbacterium,k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides; s__Bacteroides_dorei,k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides; s__Bacteroides_faecis,k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides; s__Bacteroides_fragilis,k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides; s__Bacteroides_massiliensis,k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Porphyromonadaceae; g__Parabacteroides; s__Parabacteroides_merdae,k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Bacillaceae; g__Bacillus,k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Family_XI; g__Gemella; s__Gemella_morbillorum,k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Family_XI; g__Anaerococcus; s__Anaerococcus_vaginalis,k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Family_XI; g__Ezakiella; s__Sporobacterium_sp,k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Family_XI; g__Finegoldia; s__Finegoldia_magna,k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Family_XI; g__Peptoniphilus; s__Candidatus_Peptoniphilus_massiliensis,k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Anaerostipes; s__Anaerostipes_hadrus,k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Blautia,k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Lachnoclostridium,k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae,k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Ruminococcaceae; g__Ruminococcus_1; s__Ruminococcus_bicirculans,k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Methylobacteriaceae; g__Methylobacterium,k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Cupriavidus; s__Cupriavidus_metallidurans,k__Bacteria; p__Proteobacteria; c__Deltaproteobacteria; o__Desulfovibrionales; f__Desulfovibrionaceae; g__Bilophila; s__Bilophila_wadsworthia,k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Neisseriales; f__Neisseriaceae; g__Uruburuella; s__Uruburuella_sp,k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacteriales; f__Enterobacteriaceae; g__Cronobacter; s__Cronobacter_sakazakii" \
  --o-filtered-table ../analysis/P05-clust-tree-cm-goods/$REGION/$REGION-filtered-table.qza \
  --verbose

qiime feature-table summarize \
  --i-table ../analysis/P05-clust-tree-cm-goods/$REGION/$REGION-filtered-table.qza \
  --o-visualization ../analysis/P05-clust-tree-cm-goods/$REGION/$REGION-filtered-table.qzv \
  --m-sample-metadata-file ../data/atcc-metadata.tsv \

