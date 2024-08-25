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

# Calculate the number of lines
export total_lines=$(wc -l < $INPUT_PATH/${sample_name}.marked.sam | grep -oE '[0-9]+')
echo "The total number of lines is $total_lines."
echo "========================================================================"
echo "Processing..."

# Split the file into smaller chunks with line numbers
# Create a temporary folder to store the chunks
if [ ! -d "$OUTPUT_PATH/tmp" ]; then
    mkdir -p $OUTPUT_PATH/tmp
fi

# Increase the number of chunks to improve parallelism
split -l $((total_lines / ($(nproc) * 5))) --additional-suffix=.sam $INPUT_PATH/${sample_name}.marked.sam ${OUTPUT_PATH}/tmp/${sample_name}_chunk_

# Function to process each chunk
process_chunk() {
    chunk_file=$1
    output_file=$2
    while read -r line; do
        if [[ ! $line =~ \*[\ ]*0[\ ]*0[\ ]*\*[\ ]*\*[\ ]*0[\ ]*0 ]]; then
            echo "$line" >> "$output_file"
        fi
    done < "$chunk_file"
}

export -f process_chunk

# Use GNU Parallel to process chunks in parallel with progress and increased threads
ls ${OUTPUT_PATH}/tmp/${sample_name}_chunk_*.sam | pv -l -s $(ls ${OUTPUT_PATH}/tmp/${sample_name}_chunk_*.sam | wc -l) | parallel -j $(nproc) process_chunk {} ${OUTPUT_PATH}/${sample_name}.removed.sam

# Sort the output file by line number
sort -n -k1,1 ${OUTPUT_PATH}/${sample_name}.removed.sam -o ${OUTPUT_PATH}/${sample_name}.removed.sam

echo "========================================================================"
echo "Removing unmapped complete."
# Clean up
echo "Cleaning up..."
rm -rf $OUTPUT_PATH/tmp
echo "Cleaning up complete."
