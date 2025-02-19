#! /bin/bash

# List the path value
PLINK_BASE_PATH="/mnt/raid6/bacphagenetwork/data/12_plink"
PLINK_BJ_RESULT="${PLINK_BASE_PATH}/Beijing/results/bac_age"
PLINK_GZ_RESULT="${PLINK_BASE_PATH}/Guangzhou/results/bac_age"
PLINK_FULL_RESULT="${PLINK_BASE_PATH}/Full/results/bac_age"
CURRENT_PATH=$(pwd)

echo "Initializing..."
echo "PLINK_BJ_RESULT: \t ${PLINK_BJ_RESULT}"
echo "PLINK_GZ_RESULT: \t ${PLINK_GZ_RESULT}"
echo "PLINK_FULL_RESULT: \t ${PLINK_FULL_RESULT}"
echo "CURRENT_PATH: \t ${CURRENT_PATH}"

echo "Done."
echo "========================================"

# Write all the file names to a file
echo "Writing all the file names to a file..."
echo " "
echo "Processing Beijing files..."
find "${PLINK_BJ_RESULT}" -type f | while read -r file; do
    if [[ $file == *".adjusted" ]]; then
        echo "$file" >> BJ_sbatch.list
    fi
done
echo "Done."
echo "========================================"
echo "Processing Guangzhou files..."
find "${PLINK_GZ_RESULT}" -type f | while read -r file; do
    if [[ $file == *".adjusted" ]]; then
        echo "$file" >> GZ_sbatch.list
    fi
done
echo "Done."
echo "========================================"
echo "Processing Full files again..."
find "${PLINK_FULL_RESULT}" -type f | while read -r file; do
    if [[ $file == *".adjusted" ]]; then
        echo "$file" >> FULL_sbatch.list
    fi
done
echo "Done."
echo "========================================"
echo "All files have been written to BJ_sbatch.list, GZ_sbatch.list, and FULL_sbatch.list"
