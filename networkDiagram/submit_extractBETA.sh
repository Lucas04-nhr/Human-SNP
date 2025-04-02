#! /bin/bash
#SBATCH --job-name=extractBETA
#SBATCH --output=./extractBETA.%j.out
#SBATCH --error=./extractBETA.%j.err
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

echo "Processing the progress of extracting BETA..."
export INPUT_DIR="/mnt/raid6/bacphagenetwork/data/12_plink/Full/results/bac_age"
export OUTPUT_DIR="/mnt/raid6/bacphagenetwork/data/13_networkDiagram"

echo "The input directory is $INPUT_DIR."
echo "The output directory is $OUTPUT_DIR."

# Run the Python script
python extractBETA.py