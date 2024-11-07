#! /bin/bash
#SBATCH --job-name=drawfig_BJ
#SBATCH --output=./BJ_log.%j.out
#SBATCH --error=./BJ_log.%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

# Initialize the environment
echo "Initializing the environment..."
echo "=============================="
echo ""


# Load conda
echo "Activating conda environment..."
source /home/bacphagenetwork/.bashrc
conda activate snp_analysis
echo "The conda environment has been activated."
echo "=============================="

mkdir -p /mnt/raid6/bacphagenetwork/data/12_plink/Beijing/results/fig/c005


python ../analysis.py --input /mnt/raid6/bacphagenetwork/data/12_plink/Beijing/results/c005 --output /mnt/raid6/bacphagenetwork/data/12_plink/Beijing/results/fig/c005
