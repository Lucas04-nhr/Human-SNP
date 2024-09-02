#!/bin/bash
#SBATCH --job-name=analysis_BJ
#SBATCH --output=./log/Beijing/analysis_BJ_%j.out
#SBATCH --error=./log/Beijing/analysis_BJ_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=
#SBATCH --array=1

