#!/bin/bash
#SBATCH --job-name=samtools_index_bj
#SBATCH --output=./log/samtools.%j.out
#SBATCH --error=./log/samtools.%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=ANALYSIS_PATH='/mnt/raid6/bacphagenetwork/data/bwa_analysis',INDEX_PATH='/mnt/raid6/bacphagenetwork/data/samtools_analysis',OUTPUT_PATH='/mnt/raid6/bacphagenetwork/data/samtools_result'


echo "The working directory has been changed to $ANALYSIS_PATH."

infile=($( cat bj_sbatch.list | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print $0}' ))
sample_name=$(echo "$infile" | grep -oE 'BJ[0-9]{3}')

# Convert the SAM file to BAM file
echo "Processing $sample_name..."
echo "Converting the SAM file to BAM file..."
echo "The path to the SAM file is $ANALYSIS_PATH/${sample_name}.sam."
echo "The path to the BAM file is $INDEX_PATH/${sample_name}.bam."
samtools view -bS $ANALYSIS_PATH/${sample_name}.sam > $INDEX_PATH/${sample_name}.bam || { echo "Error: samtools view failed in processing $sample_name."; exit 1; }

# Indexing the BAM file
echo "Indexing the BAM file..."
echo "The path to the BAM file is $INDEX_PATH/${sample_name}.bam."
echo "The path to the BAM index file is $OUTPUT_PATH/${sample_name}.bam.bai."
samtools index $INDEX_PATH/${sample_name}.bam $OUTPUT_PATH/${sample_name}.bam.bai || { echo "Error: samtools index failed in processing $sample_name."; exit 1; }
