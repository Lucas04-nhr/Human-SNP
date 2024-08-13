#!/bin/bash
#SBATCH --job-name=mark_dulplicates_BJ
#SBATCH --output=./log/Beijing/mark_dulplicate_BJ_%j.out
#SBATCH --error=./log/Beijing/mark_dulplicate_BJ_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=INPUT_PATH='/mnt/raid6/bacphagenetwork/data/03_samtools_sorted/Beijing',OUTPUT_PATH='/mnt/raid6/bacphagenetwork/data/04_dulplicate_marked/Beijing',JAVA_BIN='/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/bin/java',PICARD_BIN='/mnt/raid6/bacphagenetwork/tools/picard_pre-built/v3.0/picard.jar'
#SBATCH --array=1

# Initialize the environment
echo "Initializing..."

echo "The original *.bam files and their index are located in $INPUT_PATH."
echo "The marked *.bam files will be saved in $OUTPUT_PATH."

infile=($( cat bj_sbatch.list | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print $0}' ))
sample_name=$(echo "$infile" | grep -oE 'BJ[0-9]{3}')

# Check all the needed files exist

echo "Checking $INPUT_PATH/${sample_name}.sorted.bam..."
if [ ! -f "$INPUT_PATH/${sample_name}.sorted.bam" ]
then
    echo "Error: $INPUT_PATH/${sample_name}.sorted.bam does not exist."
    exit 1
else
    echo "The file $INPUT_PATH/${sample_name}.sorted.bam exists."
fi

echo "Checking $INPUT_PATH/${sample_name}.sorted.bai..."
if [ ! -f "$INPUT_PATH/${sample_name}.sorted.bai" ]
then
    echo "Error: $INPUT_PATH/${sample_name}.sorted.bai does not exist."
    exit 1
else
    echo "The file $INPUT_PATH/${sample_name}.sorted.bai exists."
fi

echo "Initializing complete."

# Mark the duplicates
echo "Marking duplicates for $INPUT_PATH/${sample_name}.sorted.bam..."

$JAVA_BIN -jar $PICARD_BIN MarkDuplicates \
    INPUT=$INPUT_PATH/${sample_name}.sorted.bam \
    OUTPUT=$OUTPUT_PATH/${sample_name}.marked.bam \
    METRICS_FILE=$OUTPUT_PATH/${sample_name}.marked.metrics.txt \
    ASSUME_SORTED=true \
    REMOVE_DUPLICATES=false \
    VALIDATION_STRINGENCY=SILENT \
|| { echo "Error: Marking duplicates for $INPUT_PATH/${sample_name}.sorted.bam failed."; exit 1; }

echo "Marking duplicates for $INPUT_PATH/${sample_name}.sorted.bam complete."
echo "The marked file is saved in $OUTPUT_PATH/${sample_name}.marked.bam."

# Index the marked file
echo "Indexing the marked file..."

samtools index -@ 4 -b $OUTPUT_PATH/${sample_name}.marked.bam $OUTPUT_PATH/${sample_name}.marked.bai || { echo "Error: Failed to index the marked file of $OUTPUT_PATH/${sample_name}.marked.bam"; exit 1; }

echo "The marked file of ${sample_name} has been indexed."
echo "The indexed file is saved in $OUTPUT_PATH/${sample_name}.marked.bai."

echo "All done."
