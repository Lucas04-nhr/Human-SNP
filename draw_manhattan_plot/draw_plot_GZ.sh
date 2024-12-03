#! /bin/bash
#SBATCH --job-name=GZ_draw
#SBATCH --output=./GZ_log.%j.out
#SBATCH --error=./GZ_log.%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

# Initialize the environment
echo "Initializing the environment..."
echo "=============================="

# Load conda
echo "Activating conda environment..."
source /home/bacphagenetwork/.bashrc
conda activate snp_analysis
echo "The conda environment has been activated."
echo "=============================="

echo "Processing the Manhattan plot of Guangzhou..."
export INPUT_DIR = "/mnt/raid6/bacphagenetwork/data/12_plink/Guangzhou/replaced/bac_age"
export OUTPUT_DIR = "/mnt/raid6/bacphagenetwork/data/12_plink/Guangzhou/merged/manhattan_plot"

echo "The input directory is /mnt/raid6/bacphagenetwork/data/12_plink/Guangzhou/replaced/bac_age."
echo "The output directory is /mnt/raid6/bacphagenetwork/data/12_plink/Guangzhou/merged/manhattan_plot."

echo "-----------------------------"

echo "Drawing the Manhattan plot..."

Rscript ./draw_plot.r --input-directory="/mnt/raid6/bacphagenetwork/data/12_plink/Guangzhou/replaced/bac_age" --output-directory="/mnt/raid6/bacphagenetwork/data/12_plink/Guangzhou/merged/manhattan_plot" \
|| { echo "Error in draw_plot.r"; exit 1; }


echo "=============================="
echo "The script has finished running."
