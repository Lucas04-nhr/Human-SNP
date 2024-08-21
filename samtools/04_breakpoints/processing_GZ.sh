#!/bin/bash
#SBATCH --job-name=mark_breakpoints_GZ
#SBATCH --output=./log/Guangzhou/mark_breakpoints_GZ_%j.out
#SBATCH --error=./log/Guangzhou/mark_breakpoints_GZ_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=INPUT_PATH='/mnt/raid6/bacphagenetwork/data/04_dulplicate_marked/Guangzhou',OUTPUT_PATH='/mnt/raid6/bacphagenetwork/data/05_breakpoint_marked/Guangzhou',JAVA_BIN='/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/bin/java',CPPFLAGS='-I/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/include',LDFLAGS='-L/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/lib/:-L/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/lib/server',GATK_OLD_BIN='/mnt/raid6/bacphagenetwork/tools/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar',GATK_NEW_BIN='/mnt/raid6/bacphagenetwork/tools/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar',PICARD_NEW_BIN='/mnt/raid6/bacphagenetwork/tools/picard_pre-built/v3.0/picard.jar',PICARD_OLD_BIN='/mnt/raid6/bacphagenetwork/tools/picard_pre-built/v2.26.0/picard.jar'
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

echo "Checking $INPUT_PATH/${sample_name}.marked.bai..."
if [ ! -f "$INPUT_PATH/${sample_name}.marked.bai" ]
then
    echo "Error: $INPUT_PATH/${sample_name}.marked.bai does not exist."
    exit 1
else
    echo "The file $INPUT_PATH/${sample_name}.marked.bai exists."
fi

echo "Initializing complete."

# Collect insert size metrics

# Mark the duplicates
echo "Marking duplicates for $INPUT_PATH/${sample_name}.sorted.bam..."

$JAVA_BIN -jar $PICARD_NEW_BIN CollectInsertSizeMetrics \
    INPUT=$INPUT_PATH/${sample_name}.marked.bam \
    OUTPUT=$OUTPUT_PATH/${sample_name}.metrics.txt \
    HISTOGRAM_FILE=$OUTPUT_PATH/${sample_name}.hist.pdf \
    MINIMUM_PCT=0.5 \
|| { echo "Error: Marking duplicates for $INPUT_PATH/${sample_name}.sorted.bam failed."; exit 1; }

echo "Marking duplicates for $INPUT_PATH/${sample_name}.marked.bam complete."
echo "The histogram is saved in $OUTPUT_PATH/${sample_name}.hist.pdf."
echo "The metrics are saved in $OUTPUT_PATH/${sample_name}.metrics.txt."

echo "All done."
