#!/bin/bash
#SBATCH --job-name=manta_exec_BJ
#SBATCH --output=./log/Beijing/manta_exec_BJ_%j.out
#SBATCH --error=./log/Beijing/manta_exec_BJ_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=PARENT_PATH='/mnt/raid6/bacphagenetwork/data/07_manta/01_exec/Beijing',MANTA_INSTALL_PATH='/mnt/raid6/bacphagenetwork/tools/manta/'
#SBATCH --array=1

# Initialize the environment
echo "Initializing..."

echo "The parent folder of the execute folders is located in $PARENT_PATH."

infile=($( cat BJ_sbatch.list | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print $0}' ))
sample_name=$(echo "$infile" | grep -oE 'BJ[0-9]{3}')

export MANTA_ANALYSIS_PATH=$PARENT_PATH/$sample_name
echo "The path to the execute folder of $sample_name has been set to $MANTA_ANALYSIS_PATH."

echo "Initializing complete."
echo "========================================"

# Run the Manta analysis
echo "Running the Manta analysis..."
${MANTA_ANALYSIS_PATH}/runWorkflow.py \
|| { echo "Error: Manta analysis failed."; exit 1; }

echo "========================================"
echo "The Manta analysis has been completed successfully."
