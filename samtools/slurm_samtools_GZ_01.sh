#!/bin/bash
#SBATCH --job-name=samtools_index_gz
#SBATCH --output=./log/Guangzhou/samtools_GZ_%j.out
#SBATCH --error=./log/Guangzhou/samtools_GZ_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=ANALYSIS_PATH='/mnt/raid6/bacphagenetwork/data/bwa_analysis/Guangzhou',INDEX_PATH='/mnt/raid6/bacphagenetwork/data/samtools_analysis/Guangzhou'
#SBATCH --array=1-160%4

# Check whether the environment exists
if conda env list | grep -q "wescall"
then
    echo "Great! The environment already exists."
    # Activate the environment
    echo "Activating the environment..."
    conda activate wescall
else
    echo "Creating the environment..."
    conda env create -n wescall -f ../requirements.txt
    echo "The environment has been created, activating it..."
    conda activate wescall
fi

echo "Initialization is complete."

echo "The working directory has been changed to $ANALYSIS_PATH."

infile=($( cat gz_01_sbatch.list | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print $0}' ))
sample_name=$(echo "$infile" | grep -oE 'GZ[0-9]{3}')

# Convert the SAM file to BAM file
echo "Processing $sample_name..."
echo "Converting the SAM file to BAM file..."
echo "The path to the SAM file is $ANALYSIS_PATH/${sample_name}.sam."
echo "The path to the BAM file is $INDEX_PATH/${sample_name}.bam."
samtools view -bS $ANALYSIS_PATH/${sample_name}.sam > $INDEX_PATH/${sample_name}.bam || { echo "Error: samtools view failed in processing $sample_name."; exit 1; }

echo "The SAM file has been successfully converted to BAM file."

# # Indexing the BAM file
# echo "Indexing the BAM file..."
# echo "The path to the BAM file is $INDEX_PATH/${sample_name}.bam."
# echo "The path to the BAM index file is $INDEX_PATH/${sample_name}.bam.bai."
# samtools index $INDEX_PATH/${sample_name}.bam $INDEX_PATH/${sample_name}.bam.bai || { echo "Error: samtools index failed in processing $sample_name."; exit 1; }

# echo "The BAM file has been successfully indexed."
