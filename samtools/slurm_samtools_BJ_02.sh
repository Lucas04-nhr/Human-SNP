#!/bin/bash
#SBATCH --job-name=samtools_sort_BJ
#SBATCH --output=./log/02/Beijing/samtools_sort_BJ_%j.out
#SBATCH --error=./log/02/Beijing/samtools_sort_BJ_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=ANALYSIS_PATH='/mnt/raid6/bacphagenetwork/data/samtools_analysis/Beijing',INDEX_PATH='/mnt/raid6/bacphagenetwork/data/samtools_index/Beijing',RESULTS_PATH='/mnt/raid6/bacphagenetwork/data/samtools_results/Beijing'
#SBATCH --array=1-201%4

conda init bash

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

# Check all the needed directories exist
echo "Checking all the needed directories exist..."
if [ ! -d "$ANALYSIS_PATH" ]
then
    mkdir -p $ANALYSIS_PATH
fi
if [ ! -d "$INDEX_PATH" ]
then
    mkdir -p $INDEX_PATH
fi
if [ ! -d "$RESULTS_PATH" ]
then
    mkdir -p $RESULTS_PATH
fi
echo "All the needed directories exist."

echo "Initialization is complete."

echo "The original *.bam files are located in $ANALYSIS_PATH."
echo "The sorted *.bam files will be saved in $RESULTS_PATH."
echo "The indexed *.bai files will be saved in $INDEX_PATH."

infile=($( cat bj_02_sbatch.list | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print $0}' ))
sample_name=$(echo "$infile" | grep -oE 'BJ[0-9]{3}')

# Check all the needed files exist
echo "Checking $ANALYSIS_PATH/${sample_name}.bam..."
if [ ! -f "$ANALYSIS_PATH/${sample_name}.bam" ]
then
    echo "Error: $ANALYSIS_PATH/${sample_name}.bam does not exist."
    exit 1
else
    echo "The file $ANALYSIS_PATH/${sample_name}.bam exists."
fi

# Sort the BAM file
echo "Processing $sample_name..."
echo "Sorting the BAM file..."
echo "The path to the BAM file is $ANALYSIS_PATH/${sample_name}.bam."
echo "The path to the sorted BAM file is $RESULTS_PATH/${sample_name}.sorted.bam."
samtools sort -@ 4 -o $RESULTS_PATH/${sample_name}.sorted.bam $ANALYSIS_PATH/${sample_name}.bam || { echo "Error: samtools sort failed in processing $sample_name."; exit 1; }

echo "The $sample_name BAM file has been successfully sorted."

# Check all the needed files exist
echo "Checking $RESULTS_PATH/${sample_name}.sorted.bam..."
if [ ! -f "$RESULTS_PATH/${sample_name}.sorted.bam" ]
then
    echo "Error: $RESULTS_PATH/${sample_name}.sorted.bam does not exist."
    exit 1
else
    echo "The file $RESULTS_PATH/${sample_name}.sorted.bam exists."
fi

# Indexing the BAM file
echo "Indexing the BAM file..."
echo "The path to the sorted BAM file is $RESULTS_PATH/${sample_name}.bam."
echo "The path to the BAM index file is $INDEX_PATH/${sample_name}.bam.bai."
samtools index $RESULTS_PATH/${sample_name}.sorted.bam $INDEX_PATH/${sample_name}.bam.bai || { echo "Error: samtools index failed in processing $sample_name."; exit 1; }

echo "The $sample_name BAM file has been successfully indexed."
