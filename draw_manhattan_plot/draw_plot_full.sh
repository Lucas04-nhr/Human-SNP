#! /bin/bash
#SBATCH --job-name=FULL_draw
#SBATCH --output=./FULL_log.%j.out
#SBATCH --error=./FULL_log.%j.err
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

echo "Processing the Manhattan plot of Beijing..."
export INPUT_DIR="/mnt/raid6/bacphagenetwork/data/12_plink/Full/modified/bac_age"
export OUTPUT_DIR="/mnt/raid6/bacphagenetwork/data/12_plink/Full/merged/manhattan_plot"

echo "The input directory is $INPUT_DIR."
echo "The output directory is $OUTPUT_DIR."

echo "-----------------------------"

echo "Drawing the Manhattan plot..."

Rscript ./draw_plot_new.r --input-directory=$INPUT_DIR --output-directory=$OUTPUT_DIR --use-saved-data=FALSE \
|| { echo "Error in draw_plot.r"; exit 1; }


echo "=============================="
echo "The script has finished running."
