#! /bin/bash

# Set the path to the genome data
GENOME_PATH='/mnt/raid6/bacphagenetwork/data/ninglab/02_samtools/03_sort'

# Check if the sbatch.list file exists, if is, ask the user if they want to overwrite it
if [ -f sbatch.list ]; then
    read -p "sbatch.list already exists. Do you want to overwrite it? (y/n) " choice
    if [[ "$choice" != "y" && "$choice" != "Y" ]]; then
      echo "Exiting without overwriting sbatch.list."
      exit 0
    else
      > sbatch.list
      echo "Overwriting sbatch.list."
    fi
fi

# Find all files ending with .bam and output their paths to a file
find "$GENOME_PATH" -type f -name '*.bam' >> sbatch.list
