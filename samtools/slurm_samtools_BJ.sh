#!/bin/bash
#SBATCH --job-name=samtools_index_bj
#SBATCH --output=./log/samtools.%j.out
#SBATCH --error=./log/samtools.%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=ANALYSIS_PATH='/mnt/raid6/bacphagenetwork/data/bwa_analysis'


echo "The working directory has been changed to $ANALYSIS_PATH."