#!/bin/bash
#SBATCH --job-name=analysis_BJ
#SBATCH --output=./log/Beijing/secondary_BJ_%j.out
#SBATCH --error=./log/Beijing/secondary_BJ_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=INPUT_PATH='/mnt/raid6/bacphagenetwork/data/06_unmapped_removed/Beijing',OUTPUT_PATH='/mnt/raid6/bacphagenetwork/data/09_secondary/Beijing',OMP_NUM_THREADS='2',DELLY_BIN
#SBATCH --array=1