#!/bin/bash
#SBATCH --job-name=sort_BJ
#SBATCH --output=./log/Beijing/sort_BJ_%j.out
#SBATCH --error=./log/Beijing/sort_BJ_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=INPUT_PATH='/mnt/raid6/bacphagenetwork/data/02_samtools_analysis/Beijing',OUTPUT_PATH='/mnt/raid6/bacphagenetwork/data/03_samtools_sorted/Beijing'
#SBATCH --array=1

# Initialize the environment
echo "Initializing..."

echo "The original *.reheader.bam files are located in $INPUT_PATH."
echo "The marked *.bam files will be saved in $OUTPUT_PATH."

infile=($( cat bj_sbatch.list | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print $0}' ))
sample_name=$(echo "$infile" | grep -oE 'BJ[0-9]{3}')

# Check all the needed files exist

echo "Checking $INPUT_PATH/${sample_name}.reheader.bam..."
if [ ! -f "$INPUT_PATH/${sample_name}.reheader.bam" ]
then
    echo "Error: $INPUT_PATH/${sample_name}.reheader.bam does not exist."
    exit 1
else
    echo "The file $INPUT_PATH/${sample_name}.reheader.bam exists."
fi

echo "Initializing complete."

# Sort the bam file
echo "Sorting the bam file of ${sample_name}..."

samtools sort -@ 4 -o $OUTPUT_PATH/${sample_name}.sorted.bam $INPUT_PATH/${sample_name}.reheader.bam || { echo "Error: Failed to sort the bam file of $INPUT_PATH/${sample_name}.bam"; exit 1; }

echo "The bam file of ${sample_name} has been sorted."
echo "The sorted file is saved in $OUTPUT_PATH/${sample_name}_sorted.bam."

# Index the sorted bam file
echo "Indexing the sorted bam file of ${sample_name}..."

samtools index -@ 4 -b $OUTPUT_PATH/${sample_name}.sorted.bam $OUTPUT_PATH/${sample_name}.sorted.bai || { echo "Error: Failed to index the sorted bam file of $OUTPUT_PATH/${sample_name}.sorted.bam"; exit 1; }

echo "The sorted bam file of ${sample_name} has been indexed."
echo "The indexed file is saved in $OUTPUT_PATH/${sample_name}.sorted.bai."

echo "All done."