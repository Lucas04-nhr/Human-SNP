#!/bin/bash
#SBATCH --job-name=bwa_analysis
#SBATCH --output=bwa_analysis.out
#SBATCH --error=bwa_analysis.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=GENOME_PATH='/mnt/raid6/bacphagenetwork/data/skin_metagenome/Beijing/02_rm_host',INDEXING_PATH='/mnt/raid6/bacphagenetwork/data/bwa_index/chm13v2.0_noY.fa',ANALYSIS_PATH='/mnt/raid6/bacphagenetwork/data/bwa_analysis'

source /etc/profile.d/modules.sh
module load anaconda3
conda init bash
source ~/.bashrc

conda activate base

srun bash ../processing.sh