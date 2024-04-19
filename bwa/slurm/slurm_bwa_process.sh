#!/bin/bash
#SBATCH --job-name=bwa_analysis
#SBATCH --output=log/bwa_analysis.out
#SBATCH --error=log/bwa_analysis.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=GENOME_PATH='/mnt/raid6/bacphagenetwork/data/skin_metagenome/Beijing/02_rm_host',INDEXING_PATH='/mnt/raid6/bacphagenetwork/data/bwa_index/chm13v2.0_noY.fa',ANALYSIS_PATH='/mnt/raid6/bacphagenetwork/data/bwa_analysis'

conda init bash
source ~/.bashrc

conda activate base

srun ../processing.sh
