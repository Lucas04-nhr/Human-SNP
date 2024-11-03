#! /bin/bash
#SBATCH --job-name=drawfig_GZ
#SBATCH --output=./GZ_log.%j.out
#SBATCH --error=./GZ_log.%j.err
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


python ../analysis.py --input /mnt/raid6/bacphagenetwork/data/12_plink/Guangzhou/results/c005/result.P1.assoc.linear --output /mnt/raid6/bacphagenetwork/data/12_plink/Guangzhou/results/fig/c005
