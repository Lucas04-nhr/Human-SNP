#!/bin/bash
#SBATCH --job-name=breakdancer_BJ
#SBATCH --output=./log/Beijing/breakdancer_BJ_%j.out
#SBATCH --error=./log/Beijing/breakdancer_BJ_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=INPUT_PATH='/mnt/raid6/bacphagenetwork/data/04_dulplicate_marked/Beijing',OUTPUT_PATH='/mnt/raid6/bacphagenetwork/data/05_breakpoint_marked/Beijing',JAVA_BIN='/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/bin/java',CPPFLAGS='-I/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/include',LDFLAGS='-L/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/lib/:-L/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/lib/server',GATK_OLD_BIN='/mnt/raid6/bacphagenetwork/tools/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar',GATK_NEW_BIN='/mnt/raid6/bacphagenetwork/tools/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar',PICARD_NEW_BIN='/mnt/raid6/bacphagenetwork/tools/picard_pre-built/v3.0/picard.jar',PICARD_OLD_BIN='/mnt/raid6/bacphagenetwork/tools/picard_pre-built/v2.26.0/picard.jar',BREAKDANCER_BIN='/mnt/raid6/bacphagenetwork/tools/breakdancer/build/bin/breakdancer-max',BREAKDANCER_CONFIG='/mnt/raid6/bacphagenetwork/tools/breakdancer/perl/bam2cfg.pl'
#SBATCH --array=1

# Initialize the environment
echo "Initializing..."

echo "The original *.bam files and their index are located in $INPUT_PATH."
echo "The marked *.bam files will be saved in $OUTPUT_PATH."

infile=($( cat bj_sbatch.list | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print $0}' ))
sample_name=$(echo "$infile" | grep -oE 'BJ[0-9]{3}')

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

# Run breakdancer

echo "Processing $INPUT_PATH/${sample_name}.marked.bam..."

# Generate the configuration file
echo "Generating the configuration file..."
$BREAKDANCER_CONFIG $INPUT_PATH/${sample_name}.marked.bam > $OUTPUT_PATH/${sample_name}.cfg \
|| { echo "Error: Generating the configuration file for $INPUT_PATH/${sample_name}.marked.bam failed."; exit 1; }

# Run breakdancer
echo "Running breakdancer..."
$BREAKDANCER_BIN $OUTPUT_PATH/${sample_name}.cfg > $OUTPUT_PATH/${sample_name}.breakdancer.txt \
|| { echo "Error: Running breakdancer for $INPUT_PATH/${sample_name}.marked.bam failed."; exit 1; }


echo "All done."
echo "The output is saved in $OUTPUT_PATH/${sample_name}.breakdancer.txt."