#!/bin/bash
#SBATCH --job-name=remove_dulplicates_BJ
#SBATCH --output=./log/Beijing/remove_dulplicates_BJ_%j.out
#SBATCH --error=./log/Beijing/remove_dulplicates_BJ_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=INPUT_PATH='/mnt/raid6/bacphagenetwork/data/05_format_converted/Beijing',OUTPUT_PATH='/mnt/raid6/bacphagenetwork/data/06_dulplicates_removed/Beijing'
#SBATCH --array=1