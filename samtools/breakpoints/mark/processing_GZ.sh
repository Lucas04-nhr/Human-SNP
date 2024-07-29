#!/bin/bash
#SBATCH --job-name=mark_dulplicates_GZ
#SBATCH --output=./log/Guangzhou/mark_dulplicates_GZ_%j.out
#SBATCH --error=./log/Guangzhou/mark_dulplicates_GZ_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=INPUT_PATH='/mnt/raid6/bacphagenetwork/data/03_samtools_sorted/Guangzhou',OUTPUT_PATH='/mnt/raid6/bacphagenetwork/data/04_samtools_marked/Guangzhou',STATS_PATH='/mnt/raid6/bacphagenetwork/data/04_samtools_marked_stats',INDEX_FILE='/mnt/raid6/bacphagenetwork/data/00_bwa_index/chm13v2/chm13v2.0_noY.fa'
#SBATCH --array=1-201%4

# Initialize the environment
echo "Initializing..."

echo "The original *.bam files are located in $INPUT_PATH."
echo "The marked *.bam files will be saved in $OUTPUT_PATH."

infile=($( cat GZ_sbatch.list | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print $0}' ))
sample_name=$(echo "$infile" | grep -oE 'GZ[0-9]{3}')

# Check all the needed files exist

echo "Checking $INPUT_PATH/${sample_name}_sorted.bam..."
if [ ! -f "$INPUT_PATH/${sample_name}_sorted.bam" ]
then
    echo "Error: $INPUT_PATH/${sample_name}_sorted.bam does not exist."
    exit 1
else
    echo "The file $INPUT_PATH/${sample_name}_sorted.bam exists."
fi

echo "Initializing complete."

# Mark the duplicates
echo "Marking the duplicates of ${sample_name}..."

samtools markdup -r -s $INPUT_PATH/${sample_name}_sorted.bam $OUTPUT_PATH/${sample_name}_marked.bam -f $STATS_PATH/${sample_name}.txt -O --reference $INDEX_FILE --write-index|| { echo "Error: Failed to mark the duplicates of $INPUT_PATH/${sample_name}.bam"; exit 1; }

echo "The duplicates of ${sample_name} have been marked."
echo "The marked file is saved in $OUTPUT_PATH/${sample_name}_marked.bam."
