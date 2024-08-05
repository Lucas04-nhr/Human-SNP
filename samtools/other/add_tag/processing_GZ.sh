#!/bin/bash
#SBATCH --job-name=add_tag_GZ
#SBATCH --output=./log/Guangzhou/add_tag_GZ_%j.out
#SBATCH --error=./log/Guangzhou/add_tag_GZ_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=INPUT_PATH='/mnt/raid6/bacphagenetwork/data/03_samtools_sorted/Guangzhou',OUTPUT_PATH='/mnt/raid6/bacphagenetwork/data/04_samtools_tag_added/Guangzhou'
#SBATCH --array=1-160%4

# Initialize the environment
echo "Initializing..."

echo "The original *.bam files are located in $INPUT_PATH."
echo "The marked *.bam files will be saved in $OUTPUT_PATH."

infile=($( cat GZ_sbatch.list | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print $0}' ))
sample_name=$(echo "$infile" | grep -oE 'GZ[0-9]{3}')

# Check all the needed files exist

echo "Checking $INPUT_PATH/${sample_name}.bam..."
if [ ! -f "$INPUT_PATH/${sample_name}_sorted.bam" ]
then
    echo "Error: $INPUT_PATH/${sample_name}_sorted.bam does not exist."
    exit 1
else
    echo "The file $INPUT_PATH/${sample_name}_sorted.bam exists."
fi

echo "Initializing complete."

# Add the tag

samtools fixmate -@ 4 -m $INPUT_PATH/${sample_name}_sorted.bam $OUTPUT_PATH/${sample_name}.bam || { echo "Error: Failed to add the tag of $INPUT_PATH/${sample_name}_sorted.bam"; exit 1; }

echo "The tag of $INPUT_PATH/${sample_name}_sorted.bam has been added successfully."
echo "The new file is saved as $OUTPUT_PATH/${sample_name}.bam."
