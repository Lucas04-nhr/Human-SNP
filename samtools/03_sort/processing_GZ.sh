#!/bin/bash
#SBATCH --job-name=sort_GZ
#SBATCH --output=./log/Guangzhou/sort_GZ_%j.out
#SBATCH --error=./log/Guangzhou/sort_GZ_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=INPUT_PATH='/mnt/raid6/bacphagenetwork/data/02_samtools_analysis/Guangzhou',OUTPUT_PATH='/mnt/raid6/bacphagenetwork/data/03_samtools_sorted/Guangzhou'
#SBATCH --array=1-160%4

# Initialize the environment
echo "Initializing..."

echo "The original *.bam files are located in $INPUT_PATH."
echo "The marked *.bam files will be saved in $OUTPUT_PATH."

infile=($( cat GZ_sbatch.list | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print $0}' ))
sample_name=$(echo "$infile" | grep -oE 'GZ[0-9]{3}')

# Check all the needed files exist

echo "Checking $INPUT_PATH/${sample_name}.bam..."
if [ ! -f "$INPUT_PATH/${sample_name}.bam" ]
then
    echo "Error: $INPUT_PATH/${sample_name}.bam does not exist."
    exit 1
else
    echo "The file $INPUT_PATH/${sample_name}.bam exists."
fi

echo "Initializing complete."

# Sort the bam file
echo "Sorting the bam file of ${sample_name}..."

samtools sort -@ 4 -o $OUTPUT_PATH/${sample_name}_sorted.bam $INPUT_PATH/${sample_name}.bam || { echo "Error: Failed to sort the bam file of $INPUT_PATH/${sample_name}.bam"; exit 1; }

echo "The bam file of ${sample_name} has been sorted."
echo "The sorted file is saved in $OUTPUT_PATH/${sample_name}_sorted.bam."

echo "All done."