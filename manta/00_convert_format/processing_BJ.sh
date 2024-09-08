#!/bin/bash
#SBATCH --job-name=analysis_BJ
#SBATCH --output=./log/Beijing/converting_format_BJ_%j.out
#SBATCH --error=./log/Beijing/converting_format_BJ_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=INPUT_PATH='/mnt/raid6/bacphagenetwork/data/06_unmapped_removed/Beijing',OUTPUT_PATH='/mnt/raid6/bacphagenetwork/data/07_manta/00_format_converted/Beijing'
#SBATCH --array=2-201%5

# Initialize the environment
echo "Initializing..."

echo "The original *.sam files are located in $INPUT_PATH."
echo "The converted *.bam files and their indexes will be saved in $OUTPUT_PATH."

infile=($( cat BJ_sbatch.list | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print $0}' ))
sample_name=$(echo "$infile" | grep -oE 'BJ[0-9]{3}')

# Check all the needed files exist
# Check the output folder
if [ ! -d "$OUTPUT_PATH" ]
then
    echo "The output folder does not exist."
    exit 1
fi

# Check the input file
echo "Checking $INPUT_PATH/${sample_name}.removed.sam..."
if [ ! -f "$INPUT_PATH/${sample_name}.removed.sam" ]
then
    echo "Error: $INPUT_PATH/${sample_name}.removed.sam does not exist."
    exit 2
else
    echo "The file $INPUT_PATH/${sample_name}.removed.sam exists."
fi

# Check the output file
echo "Checking $OUTPUT_PATH/${sample_name}.bam..."
if [ -f "$OUTPUT_PATH/${sample_name}.bam" ]
then
    echo "Warning: $OUTPUT_PATH/${sample_name}.removed.sam exists and will be overwritten."
    rm $OUTPUT_PATH/${sample_name}.bam
fi

echo "Checking $OUTPUT_PATH/${sample_name}.bam.bai..."
if [ -f "$OUTPUT_PATH/${sample_name}.bam.bai" ]
then
    echo "Warning: $OUTPUT_PATH/${sample_name}.removed.sam exists and will be overwritten."
    rm $OUTPUT_PATH/${sample_name}.bam.bai
fi

echo "Initializing complete."
echo "=========================================================================="
echo "Processing ${sample_name}..."

# Convert the sam file to bam file
echo "Converting the sam file to bam file..."

samtools view -@ 4 -bS $INPUT_PATH/${sample_name}.removed.sam > $OUTPUT_PATH/${sample_name}.bam \
|| { echo "Error: samtools view failed"; exit 3; }

echo "Converting the sam file to bam file complete."

# Index the bam file
echo "Indexing the bam file..."

samtools index -@ 4 $OUTPUT_PATH/${sample_name}.bam \
|| { echo "Error: samtools index failed"; exit 4; }

echo "Indexing the bam file complete."

echo "Processing ${sample_name} complete."
echo "=========================================================================="

echo "The job has been completed successfully."
echo "The output files are located in $OUTPUT_PATH."
