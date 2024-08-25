#!/bin/bash
#SBATCH --job-name=remove_unmapped_BJ
#SBATCH --output=./log/Beijing/remove_unmapped_BJ_%j.out
#SBATCH --error=./log/Beijing/remove_unmapped_BJ_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=INPUT_PATH='/mnt/raid6/bacphagenetwork/data/05_format_converted/Beijing',OUTPUT_PATH='/mnt/raid6/bacphagenetwork/data/06_unmapped_removed/Beijing'
#SBATCH --array=1


# Initialize the environment
echo "Initializing..."

echo "The original *.sam files are located in $INPUT_PATH."
echo "The marked *.sam files will be saved in $OUTPUT_PATH."

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
echo "Checking $INPUT_PATH/${sample_name}.marked.sam..."
if [ ! -f "$INPUT_PATH/${sample_name}.marked.sam" ]
then
    echo "Error: $INPUT_PATH/${sample_name}.marked.sam does not exist."
    exit 1
else
    echo "The file $INPUT_PATH/${sample_name}.marked.sam exists."
fi

echo "Initializing complete."

# Remove unmapped
echo "Removing unmapped from $INPUT_PATH/${sample_name}.marked.sam to $OUTPUT_PATH/${sample_name}.removed.sam..."

i = 0

while read -r line; do
    i = $((i+1))
    echo "Processing line $i..."
    if [[ ! $line =~ \*[\ ]*0[\ ]*0[\ ]*\*[\ ]*\*[\ ]*0[\ ]*0 ]]; then
        echo "$line" >> $OUTPUT_PATH/${sample_name}.removed.sam
    fi
    echo "Line $i processed."
done < $INPUT_PATH/${sample_name}.marked.sam

echo "Removing unmapped complete."
