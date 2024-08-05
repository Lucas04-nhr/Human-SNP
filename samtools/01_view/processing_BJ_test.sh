#!/bin/bash
#SBATCH --job-name=view_and_index_BJ
#SBATCH --output=./log/Beijing/view_and_index_BJ_%j.out
#SBATCH --error=./log/Beijing/view_and_index_BJ_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=IMPUT_PATH='/mnt/raid6/bacphagenetwork/data/01_bwa_analysis/Beijing',OUTPUT_PATH='/mnt/raid6/bacphagenetwork/data/02_samtools_viewed/Beijing'
#SBATCH --array=1

echo "Initialization is complete."

echo "The working directory has been changed to $IMPUT_PATH."

infile=($( cat BJ_sbatch.list | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print $0}' ))
sample_name=$(echo "$infile" | grep -oE 'BJ[0-9]{3}')

# Convert the SAM file to BAM file
echo "Processing $sample_name..."
echo "Converting the SAM file to BAM file..."
echo "The path to the SAM file is $IMPUT_PATH/${sample_name}.sam."
echo "The path to the BAM file is $OUTPUT_PATH/${sample_name}.bam."
samtools view -bHS $IMPUT_PATH/${sample_name}.sam > $OUTPUT_PATH/${sample_name}.bam || { echo "Error: samtools view failed in processing $sample_name."; exit 1; }

echo "The SAM file $sample_name has been successfully converted to BAM file."
