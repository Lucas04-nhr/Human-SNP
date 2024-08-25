#!/bin/bash
#SBATCH --job-name=remove_unmapped_GZ
#SBATCH --output=./log/Guangzhou/remove_unmapped_GZ_%j.out
#SBATCH --error=./log/Guangzhou/remove_unmapped_GZ_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=INPUT_PATH='/mnt/raid6/bacphagenetwork/data/05_format_converted/Guangzhou',OUTPUT_PATH='/mnt/raid6/bacphagenetwork/data/06_unmapped_removed/Guangzhou'
#SBATCH --array=2-201%5


# Initialize the environment
echo "Initializing..."

echo "The original *.sam files are located in $INPUT_PATH."
echo "The marked *.sam files will be saved in $OUTPUT_PATH."

infile=($( cat GZ_sbatch.list | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print $0}' ))
sample_name=$(echo "$infile" | grep -oE 'GZ[0-9]{3}')

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

# Check the output file
echo "Checking $OUTPUT_PATH/${sample_name}.removed.sam..."
if [ -f "$OUTPUT_PATH/${sample_name}.removed.sam" ]
then
    echo "Warning: $OUTPUT_PATH/${sample_name}.removed.sam exists and will be overwritten."
    rm $OUTPUT_PATH/${sample_name}.removed.sam
fi

echo "Initializing complete."

# Remove unmapped
echo "Removing unmapped from $INPUT_PATH/${sample_name}.marked.sam to $OUTPUT_PATH/${sample_name}.removed.sam..."

# Using samtools
echo "Processing ${sample_name}..."
samtools view -F 4 -h $INPUT_PATH/${sample_name}.marked.sam > $OUTPUT_PATH/${sample_name}.removed.sam \
|| { echo "Error: Processing failed for ${sample_name}."; exit 1; }
echo "Processing complete for ${sample_name}."


# # Calculate the number of lines
# export total_lines=$(wc -l < $INPUT_PATH/${sample_name}.marked.sam | grep -oE '[0-9]+')
# echo "The total number of lines is $total_lines."
# echo "========================================================================"
# echo "Processing..."

# export i=0

# while read -r line; do
#     export i=$((i+1))
#     echo "Processing line $i of $total_lines..."
#     if [[ ! $line =~ \*[\ ]*0[\ ]*0[\ ]*\*[\ ]*\*[\ ]*0[\ ]*0 ]]; then
#         echo "$line" >> $OUTPUT_PATH/${sample_name}.removed.sam
#         echo "Line $i will be saved."
#     else
#         echo "Line $i will be removed."
#     fi
# done < $INPUT_PATH/${sample_name}.marked.sam

# echo "========================================================================"
# echo "Removing unmapped complete."
