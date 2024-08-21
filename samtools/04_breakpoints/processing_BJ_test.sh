#!/bin/bash
#SBATCH --job-name=mark_breakpoints_BJ
#SBATCH --output=./log/Beijing/mark_breakpoints_BJ_%j.out
#SBATCH --error=./log/Beijing/mark_breakpoints_BJ_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=INPUT_PATH='/mnt/raid6/bacphagenetwork/data/04_dulplicate_marked/Beijing',OUTPUT_PATH='/mnt/raid6/bacphagenetwork/data/05_breakpoint_marked/Beijing',JAVA_BIN='/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/bin/java',PATH='/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/bin/java:$PATH',CPPFLAGS='-I/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/include',LDFLAGS='-L/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/lib/:-L/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/lib/server',GATK_OLD_BIN='/mnt/raid6/bacphagenetwork/tools/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar',GATK_NEW_BIN='/mnt/raid6/bacphagenetwork/tools/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar'
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



echo "All done."
