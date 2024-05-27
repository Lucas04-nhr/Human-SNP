#!/bin/bash
#SBATCH --job-name=samtools_sort_GZ
#SBATCH --output=./log/Guangzhou/samtools_sort_GZ_%j.out
#SBATCH --error=./log/Guangzhou/samtools_sort_GZ_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=ANALYSIS_PATH='/mnt/raid6/bacphagenetwork/data/samtools_analysis/Guangzhou',INDEX_PATH='/mnt/raid6/bacphagenetwork/data/samtools_index/Guangzhou',RESULTS_PATH='/mnt/raid6/bacphagenetwork/data/samtools_results/Guangzhou'
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

echo "The original *.bam files are located in $ANALYSIS_PATH."
echo "The sorted *.bam files will be saved in $RESULTS_PATH."
echo "The indexed *.bai files will be saved in $INDEX_PATH."

infile=($( cat gz_sbatch.list | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print $0}' ))
sample_name=$(echo "$infile" | grep -oE 'GZ[0-9]{3}')

# Sort the BAM file
echo "Processing $sample_name..."
echo "Sorting the BAM file..."
echo "The path to the BAM file is $ANALYSIS_PATH/${sample_name}.bam."
echo "The path to the sorted BAM file is $RESULTS_PATH/${sample_name}.sorted.bam."
samtools sort -@ 4 -o $RESULTS_PATH/${sample_name}.sorted.bam $ANALYSIS_PATH/${sample_name}.bam || { echo "Error: samtools sort failed in processing $sample_name."; exit 1; }


# Indexing the BAM file
echo "Indexing the BAM file..."
echo "The path to the sorted BAM file is $RESULTS_PATH/${sample_name}.bam."
echo "The path to the BAM index file is $INDEX_PATH/${sample_name}.bam.bai."
samtools index $RESULTS_PATH/${sample_name}.bam $INDEX_PATH/${sample_name}.bam.bai || { echo "Error: samtools index failed in processing $sample_name."; exit 1; }

echo "The $sample_name BAM file has been successfully indexed."
