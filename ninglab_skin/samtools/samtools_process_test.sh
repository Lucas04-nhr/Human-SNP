#! /bin/bash
#SBATCH --job-name=ninglab_data_samtools
#SBATCH --output=./log/ninglab_data_samtools.%j.out
#SBATCH --error=./log/ninglab_data_samtools.%j.err
#SBATCH --cpus-per-task=4
#SBATCH --array=1
#SBATCH --mem=32G

# Initialize the environment
echo "Initializing the environment..."
echo "=============================="

# Load conda
echo "Activating conda environment..."
source /home/bacphagenetwork/.bashrc
conda activate snp_analysis
echo "The conda environment has been activated."
echo "=============================="


# Set the path to the genome data
GENOME_PATH='/mnt/raid6/bacphagenetwork/data/ninglab_skin'
BASE_PATH='/mnt/raid6/bacphagenetwork/data/ninglab'
INDEXING_PATH='/mnt/raid6/bacphagenetwork/data/00_bwa_index/GRCh38'
INDEXING_FILE='/mnt/raid6/bacphagenetwork/data/00_bwa_index/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa'
ANALYSIS_PATH="${BASE_PATH}/01_bwa_analysis"
SAMTOOLS_PATH="${BASE_PATH}/02_samtools"

# Prompt the paths
echo "The path to the genome data has been set to $GENOME_PATH."
echo "The indexing file is located at $INDEXING_FILE."
echo "The analysis path is set to $ANALYSIS_PATH."
echo "The samtools path is set to $SAMTOOLS_PATH."
echo "Checking whether the samtools result file path exists..."
  
if [ -d "$SAMTOOLS_PATH" ]; then
    echo "The samtools result file path exists."
else
    echo "The samtools result file path does not exist. Creating"
    mkdir -p "$SAMTOOLS_PATH/01_view"
    mkdir -p "$SAMTOOLS_PATH/02_rehead"
    mkdir -p "$SAMTOOLS_PATH/03_sort"
    echo "The samtools result file path has been created."
fi
echo "=============================="

# Analysing the data
echo "Using samtools to analyse the data..."
infile=($( cat sbatch.list | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print $0}' ))
sample_name=$(basename "${infile}" .sam)

echo "Processing ${sample_name}..."
echo "The path to the SAM file is $ANALYSIS_PATH/${sample_name}.sam."
# Check if the SAM file exists
if [ ! -f "$ANALYSIS_PATH/${sample_name}.sam" ]; then
    echo "Error: The SAM file $ANALYSIS_PATH/${sample_name}.sam does not exist."
    exit 1
else
    echo "The SAM file $ANALYSIS_PATH/${sample_name}.sam exists."
fi
echo "==============================="

# Convert the SAM file to BAM file
echo "Converting the SAM file to BAM file..."
samtools view -bS "$ANALYSIS_PATH/${sample_name}.sam" > "$SAMTOOLS_PATH/01_view/${sample_name}.bam" \
|| { echo "Error: SAM to BAM conversion failed."; exit 1; }
echo "The SAM file ${sample_name} has been successfully converted to BAM file."
echo "================================"

# Create the header file
echo "Creating the header file..."
samtools view -H "$SAMTOOLS_PATH/01_view/${sample_name}.bam" > "$SAMTOOLS_PATH/01_view/${sample_name}.header" \
|| { echo "Error: Header file creation failed."; exit 1; }
echo "The header file has been created."

# Adding some lines to the header file
echo "Adding some lines to the header file..."
echo -e "@RG\tID:${sample_name}\tSM:${sample_name}\tPL:ILLUMINA" >> "$SAMTOOLS_PATH/01_view/${sample_name}.header" \
|| { echo "Error: Failed to add lines to the header file."; exit 1; }
echo -e "@HD\tVN:1.0\tSO:coordinate" >> "$SAMTOOLS_PATH/01_view/${sample_name}.header" \
|| { echo "Error: Failed to add lines to the header file."; exit 1; }
echo "Some lines have been added to the header file."

# Combine the header and the BAM file
echo "Reheadering the BAM file..."
samtools reheader "$SAMTOOLS_PATH/01_view/${sample_name}.header" "$SAMTOOLS_PATH/01_view/${sample_name}.bam" > "$SAMTOOLS_PATH/02_rehead/${sample_name}.reheader.bam" \
|| { echo "Error: BAM reheadering failed."; exit 1; }
echo "The BAM file has been reheadered."
echo "The reheadered BAM file is saved in $SAMTOOLS_PATH/02_rehead/${sample_name}.reheader.bam."
echo "==============================="

# Sort the BAM file
echo "Sorting the BAM file..."
samtools sort -@ 4 -o "$SAMTOOLS_PATH/03_sort/${sample_name}.sorted.bam" "$SAMTOOLS_PATH/02_rehead/${sample_name}.reheader.bam" \
|| { echo "Error: BAM sorting failed."; exit 1; }
echo "The BAM file has been sorted."
echo "The sorted BAM file is saved in $SAMTOOLS_PATH/03_sort/${sample_name}.sorted.bam."

# Index the sorted BAM file
echo "Indexing the sorted BAM file..."
samtools index -@ 4 -b "$SAMTOOLS_PATH/03_sort/${sample_name}.sorted.bam" "$SAMTOOLS_PATH/03_sort/${sample_name}.sorted.bai" \
|| { echo "Error: BAM indexing failed."; exit 1; }
echo "The sorted BAM file has been indexed."
echo "================================"
echo "All done."
echo "The analysis of ${sample_name} has been completed."
echo "Here are the output files:"
echo "1. $SAMTOOLS_PATH/01_view/${sample_name}.bam"
echo "2. $SAMTOOLS_PATH/01_view/${sample_name}.header"
echo "3. $SAMTOOLS_PATH/02_rehead/${sample_name}.reheader.bam"
echo "4. $SAMTOOLS_PATH/03_sort/${sample_name}.sorted.bam"
echo "5. $SAMTOOLS_PATH/03_sort/${sample_name}.sorted.bai"
