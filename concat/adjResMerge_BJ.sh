#! /bin/bash
#SBATCH --job-name=BJ_merge
#SBATCH --output=./log/BJ_merge_log.%j.out
#SBATCH --error=./log/BJ_merge_log.%j.err
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

# Merge the files
echo "Merging the files of Beijing..."
python3 ./adjResMerge.py --input-directory=/mnt/raid6/bacphagenetwork/data/12_plink/Beijing/results/modified/c004 \
|| { echo "Error during the merging process."; exit 1; }
echo "The files have been merged."

echo "=============================="
echo "The script has finished running."
