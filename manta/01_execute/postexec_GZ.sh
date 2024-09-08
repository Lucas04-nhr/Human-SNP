#!/bin/bash
#SBATCH --job-name=manta_copy_result_GZ
#SBATCH --output=./log/Guangzhou/02/manta_copy_result_GZ_%j.out
#SBATCH --error=./log/Guangzhou/02/manta_copy_result_GZ_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=PARENT_PATH='/mnt/raid6/bacphagenetwork/data/07_manta/01_exec/Guangzhou',MANTA_INSTALL_PATH='/mnt/raid6/bacphagenetwork/tools/manta',RESULT_PATH='/mnt/raid6/bacphagenetwork/niehaoran/Human-SNP/manta/02_result/Guangzhou'
#SBATCH --array=1-160%5
#SBATCH --dependency=afterok:59100

# Initialize the environment
echo "Initializing..."

echo "The parent folder of the execute folders is located in $PARENT_PATH."

infile=($( cat GZ_sbatch.list | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print $0}' ))
sample_name=$(echo "$infile" | grep -oE 'GZ[0-9]{3}')

export MANTA_ANALYSIS_PATH=$PARENT_PATH/$sample_name
echo "The path to the execute folder of $sample_name has been set to $MANTA_ANALYSIS_PATH."

echo "Initializing complete."
echo "========================================"

# Copy the Manta analysis result
echo "Copying the Manta analysis result of ${sample_name}..."
export STATS_PATH=${MANTA_ANALYSIS_PATH}/results/stats

if [ ! -d "${RESULT_PATH}/${sample_name}" ]
then
    mkdir -p ${RESULT_PATH}/${sample_name}
    echo "The result folder has been created."
else
    echo "The result folder already exists,overwriting it..."
    rm -rf ${RESULT_PATH}/${sample_name}
fi

cp -r ${STATS_PATH} ${RESULT_PATH} \
|| { echo "Error: Copying the Manta analysis result failed."; exit 1; }
echo "Copying complete."

echo "Renaming the result folder..."
mv ${RESULT_PATH}/stats ${RESULT_PATH}/${sample_name} \
|| { echo "Error: Renaming the result folder failed."; exit 1; }
echo "Renaming complete."

echo "========================================"
echo "The Manta analysis result has been copied successfully."
