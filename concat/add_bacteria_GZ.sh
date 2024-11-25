#! /bin/bash
#SBATCH --job-name=GZ_add_bacteria
#SBATCH --output=./log/Guangzhou/GZ_log.%j.out
#SBATCH --error=./log/Guangzhou/GZ_log.%j.err
#SBATCH --cpus-per-task=2
#SBATCH --array=1-72%4
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

infile=($( cat GZ_sbatch.list | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print $0}' ))

# Extract the Phenotype Col Number
pheno_col=$(echo $infile | grep -o -E '[0-9]+' | head -n 2)
echo "The phenotype column number is $pheno_col."
python3 ./add_bacteria.py --input-file=${infile} \
|| { echo "Error in add_bacteria.py"; exit 1; }
echo "The bacteria have been added to the file."
echo "=============================="
echo "The script has finished running."
