#!/bin/bash
#SBATCH --job-name=analysis_BJ
#SBATCH --output=./log/Beijing/secondary_BJ_%j.out
#SBATCH --error=./log/Beijing/secondary_BJ_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=INPUT_PATH='/mnt/raid6/bacphagenetwork/data/06_unmapped_removed/Beijing',OUTPUT_PATH='/mnt/raid6/bacphagenetwork/data/09_secondary/Beijing',OMP_NUM_THREADS='2',DELLY_BIN
#SBATCH --array=1

# Initialate conda
source ~/.bashrc

# conda activate snp_analysis

echo "The *.sam files whose unmapped area has been removed are located in $INPUT_PATH."
echo "The analysis files will be saved in $OUTPUT_PATH."

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
    exit 2
else
    echo "The file $INPUT_PATH/${sample_name}.removed.sam exists."
fi

# Process the file
echo "Processing ${sample_name}..."

# Convert the sam file to bam file
echo "Converting the sam file to bam file..."
if [ -f "$OUTPUT_PATH/${sample_name}.removed.bam" ]
then
    echo "The file $OUTPUT_PATH/${sample_name}.removed.bam already exists."
else
    samtools view -bS $INPUT_PATH/${sample_name}.removed.sam > $OUTPUT_PATH/${sample_name}.removed.bam \
    || { echo "Error: samtools view failed"; exit 3; }
    echo "The file $OUTPUT_PATH/${sample_name}.removed.bam has been created."
fi

# Export the config file
echo "Exporting the config file..."
echo "$OUTPUT_PATH/${sample_name}.removed.bam" > $OUTPUT_PATH/config.txt

# Run breakdancer-max
echo "Running breakdancer-max..."
$BREAKDANCER_BIN -q 10 -r 0.1 -h 200 -c $OUTPUT_PATH/config.txt > $OUTPUT_PATH/${sample_name}.bdout.txt \
|| { echo "Error: breakdancer-max failed"; exit 4; }

# Remove the intermediate files
echo "Removing the intermediate files..."
rm $OUTPUT_PATH/${sample_name}.removed.bam
echo "The intermediate files have been removed."

echo "The file $OUTPUT_PATH/${sample_name}.bdout.txt has been created."
echo "The analysis of ${sample_name} has been completed."
