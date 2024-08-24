#!/bin/bash
#SBATCH --job-name=convert_format_GZ
#SBATCH --output=./log/Guangzhou/convert_format_GZ_%j.out
#SBATCH --error=./log/Guangzhou/convert_format_GZ_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=INPUT_PATH='/mnt/raid6/bacphagenetwork/data/04_dulplicate_marked/Guangzhou',OUTPUT_PATH='/mnt/raid6/bacphagenetwork/data/05_format_converted/Guangzhou'
#SBATCH --array=1-160%5

# Initialize the environment
echo "Initializing..."

echo "The original *.bam files and their index are located in $INPUT_PATH."
echo "The marked *.bam files will be saved in $OUTPUT_PATH."

infile=($( cat gz_sbatch.list | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print $0}' ))
sample_name=$(echo "$infile" | grep -oE 'GZ[0-9]{3}')

# Check all the needed files exist

echo "Checking $INPUT_PATH/${sample_name}.marked.bam..."
if [ ! -f "$INPUT_PATH/${sample_name}.marked.bam" ]
then
    echo "Error: $INPUT_PATH/${sample_name}.marked.bam does not exist."
    exit 1
else
    echo "The file $INPUT_PATH/${sample_name}.marked.bam exists."
fi

echo "Initializing complete."

# Convert the format

echo "Converting the format of $INPUT_PATH/${sample_name}.marked.bam to $OUTPUT_PATH/${sample_name}.marked.sam..."

samtools view -@ 4 -h $INPUT_PATH/${sample_name}.marked.bam > $OUTPUT_PATH/${sample_name}.marked.sam \
|| { echo "Error: samtools view failed."; exit 1; }

echo "All done."
