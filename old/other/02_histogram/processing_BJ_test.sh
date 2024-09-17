#!/bin/bash
#SBATCH --job-name=histogram_BJ
#SBATCH --output=./log/Beijing/histogram_BJ_%j.out
#SBATCH --error=./log/Beijing/histogram_BJ_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=INPUT_PATH='/mnt/raid6/bacphagenetwork/data/06_unmapped_removed/Beijing',OUTPUT_PATH='/mnt/raid6/bacphagenetwork/data/07_histogram/Beijing',JAVA_HOME='/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/',JAVA_BIN='/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/bin/java',PICARD_BIN='/mnt/raid6/bacphagenetwork/tools/picard_pre-built/v3.0/picard.jar',PICARD_OLD_BIN='/mnt/raid6/bacphagenetwork/tools/picard_pre-built/v2.26.0/picard.jar',LDFLAGS='-L/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/lib/server',CPPFLAGS='-I/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/include'
#SBATCH --array=1


# Initialize the environment
echo "Initializing..."

echo "The *.sam files whose unmapped area has been removed are located in $INPUT_PATH."
echo "The histogram files will be saved in $OUTPUT_PATH."

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
    exit 1
else
    echo "The file $INPUT_PATH/${sample_name}.removed.sam exists."
fi

# Draw the histogram

echo "Drawing the histogram for $INPUT_PATH/${sample_name}.removed.sam..."

$JAVA_BIN -jar $PICARD_BIN CollectInsertSizeMetrics \
    I=$INPUT_PATH/${sample_name}.removed.sam \
    O=$OUTPUT_PATH/${sample_name}.histogram.txt \
    H=$OUTPUT_PATH/${sample_name}.histogram.pdf \
    M=0.5 \
|| { echo "Error: Drawing the histogram for $INPUT_PATH/${sample_name}.removed.sam failed."; exit 1; }

echo "Drawing the histogram for $INPUT_PATH/${sample_name}.removed.sam complete."
echo "The histogram is saved in $OUTPUT_PATH/${sample_name}.histogram.pdf."
echo "The histogram data is saved in $OUTPUT_PATH/${sample_name}.histogram.txt."

echo "All done."