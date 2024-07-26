#!/bin/bash
#SBATCH --job-name=samtools_collate_BJ
#SBATCH --output=./log/Beijing/samtools_sort_BJ_%j.out
#SBATCH --error=./log/Beijing/samtools_sort_BJ_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=INPUT_PATH='/mnt/raid6/bacphagenetwork/data/02_samtools_analysis/Beijing',OUTPUT_PATH='/mnt/raid6/bacphagenetwork/data/03_samtools_collate/Beijing'
#SBATCH --array=1-201%4

echo "The original *.bam files are located in $INPUT_PATH."
echo "The collated *.bed files will be saved in $OUTPUT_PATH."

infile=($( cat bj_sbatch.list | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print $0}' ))
sample_name=$(echo "$infile" | grep -oE 'BJ[0-9]{3}')

# Check all the needed files exist
echo "Checking $INPUT_PATH/${sample_name}.bam..."
if [ ! -f "$INPUT_PATH/${sample_name}.bam" ]
then
    echo "Error: $INPUT_PATH/${sample_name}.bam does not exist."
    exit 1
else
    echo "The file $INPUT_PATH/${sample_name}.bam exists."
fi

# Perform the collation
echo "Performing the collation..."
samtools collate -u $INPUT_PATH/${sample_name}.bam | awk '{if ($9 > 0 && $9 < 1000) print $1, $3, $4, $9}' > $OUTPUT_PATH/${sample_name}.bed || { echo "Error: Failed to collate $INPUT_PATH/${sample_name}.bam"; exit 1; }

echo "The collation of $INPUT_PATH/${sample_name}.bam has been completed."