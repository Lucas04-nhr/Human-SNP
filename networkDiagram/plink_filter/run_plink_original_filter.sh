#! /bin/bash
#SBATCH --job-name=plink_original_filter
#SBATCH --output=plink_original_filter.%j.out
#SBATCH --error=plink_original_filter.%j.err
#SBATCH --cpus-per-task=1
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

echo "Running plink_original_filter.py script..."
export INPUT_DIR="/mnt/raid6/bacphagenetwork/data/12_plink/Full/results/bac_age"
export OUTPUT_FILE="/mnt/raid6/bacphagenetwork/data/13_networkDiagram/plink_original_filtered.csv"

echo "The input directory is $INPUT_DIR."
echo "The output file is $OUTPUT_FILE."

# Run the Python script
python /mnt/raid6/bacphagenetwork/niehaoran/Human-SNP/networkDiagram/plink_original_filter.py -d $INPUT_DIR -o $OUTPUT_FILE