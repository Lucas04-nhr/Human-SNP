#!/bin/bash
#SBATCH --job-name=bwa_analysis
#SBATCH --output=bwa_analysis.out
#SBATCH --error=bwa_analysis.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --partition=standard

source /etc/profile.d/modules.sh
module load anaconda3
conda init bash
source ~/.bashrc

conda activate base

srun bash processing.sh