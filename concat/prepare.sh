#! /bin/bash

# List the path value
PLINK_BASE_PATH="/mnt/raid6/bacphagenetwork/data/12_plink"
PLINK_BJ_RESULT="${PLINK_BASE_PATH}/Beijing/results/original/c004"
PLINK_GZ_RESULT="${PLINK_BASE_PATH}/Guangzhou/results/original/c004"
CURRENT_PATH=$(pwd)

echo "Initializing..."
echo "PLINK_BJ_RESULT: \t ${PLINK_BJ_RESULT}"
echo "PLINK_GZ_RESULT: \t ${PLINK_GZ_RESULT}"
echo "CURRENT_PATH: \t ${CURRENT_PATH}"

echo "Done."
echo "========================================"

# Write all the file names to a file
echo "Writing all the file names to a file..."
echo " "
echo "Processing Beijing files..."
for file in ${PLINK_BJ_RESULT}; do
    echo ${file} >> ${CURRENT_PATH}/BJ_sbatch.list \
    || { echo "Error: Failed to write to BJ_sbatch.list"; exit 1; }
done
echo "Done."
echo "----------------------------------------"
echo "Processing Guangzhou files..."
for file in ${PLINK_GZ_RESULT}; do
    echo ${file} >> ${CURRENT_PATH}/GZ_sbatch.list \
    || { echo "Error: Failed to write to GZ_sbatch.list"; exit 1; }
done
echo "Done."
echo "========================================"
echo "All files have been written to BJ_sbatch.list and GZ_sbatch.list"
