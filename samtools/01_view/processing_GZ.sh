#!/bin/bash
#SBATCH --job-name=view_and_index_GZ
#SBATCH --output=./log/Guangzhou/view_and_index_GZ_%j.out
#SBATCH --error=./log/Guangzhou/view_and_index_GZ_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=IMPUT_PATH='/mnt/raid6/bacphagenetwork/data/01_bwa_analysis/Guangzhou',OUTPUT_PATH='/mnt/raid6/bacphagenetwork/data/02_samtools_viewed/Guangzhou'
#SBATCH --array=1-160%4

echo "Initialization is complete."

echo "The working directory has been changed to $IMPUT_PATH."

infile=($( cat GZ_sbatch.list | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print $0}' ))
sample_name=$(echo "$infile" | grep -oE 'GZ[0-9]{3}')

echo "Processing $sample_name..."
echo "The path to the SAM file is $IMPUT_PATH/${sample_name}.sam."
echo "The path to the BAM file is $OUTPUT_PATH/${sample_name}.bam."

# Convert the SAM file to BAM file
echo "Converting the SAM file to BAM file..."
samtools view -bS $IMPUT_PATH/${sample_name}.sam > $OUTPUT_PATH/${sample_name}.bam || { echo "Error: samtools view failed in processing $sample_name."; exit 1; }
echo "The SAM file $sample_name has been successfully converted to BAM file."

# Create the header file
echo "Creating the header file..."
samtools view -H $OUTPUT_PATH/${sample_name}.bam > $OUTPUT_PATH/${sample_name}.header || { echo "Error: samtools view failed in creating the header file for $sample_name."; exit 1; }
echo "The header file has been created."

# Adding some lines to the header file
echo "Adding some lines to the header file..."
echo -e "@RG\tID:${sample_name}\tSM:${sample_name}\tPL:ILLUMINA" >> $OUTPUT_PATH/${sample_name}.header || { echo "Error: Failed to add lines to the header file for $sample_name."; exit 1; }
echo -e "@HD\tVN:1.0\tSO:coordinate" >> $OUTPUT_PATH/${sample_name}.header || { echo "Error: Failed to add lines to the header file for $sample_name."; exit 1; }
echo "Some lines have been added to the header file."

# Combine the header and the BAM file
echo "Reheadering the BAM file..."
samtools reheader $OUTPUT_PATH/${sample_name}.header $OUTPUT_PATH/${sample_name}.bam > $OUTPUT_PATH/${sample_name}.reheader.bam || { echo "Error: samtools reheader failed in reheadering the BAM file for $sample_name."; exit 1; }
echo "The BAM file has been reheadered."
echo "The processing of $sample_name has been completed."
echo "The output file is $OUTPUT_PATH/${sample_name}.reheader.bam."
