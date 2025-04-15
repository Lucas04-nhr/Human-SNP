#! /bin/bash
#SBATCH --job-name=ninglab_data_gatk
#SBATCH --output=./log/ninglab_data_gatk.%j.out
#SBATCH --error=./log/ninglab_data_gatk.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lucas04@hust.edu.cn
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
INDEXING_FILE='/mnt/raid6/bacphagenetwork/data/00_bwa_index/GRCh38/.original/Homo_sapiens.GRCh38.dna.toplevel.fa'
ANALYSIS_PATH="${BASE_PATH}/01_bwa_analysis"
SAMTOOLS_PATH="${BASE_PATH}/02_samtools"
GATK_PATH="${BASE_PATH}/03_gatk"

# Set other environment variables
KNOWN_SITES_PATH="$BASE_PATH/../00_bwa_index/GRCh38/known-sites/dbsnp138"
KNOWN_SITES_FILE="$KNOWN_SITES_PATH/hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"
SORTED_DATA_PATH="$SAMTOOLS_PATH/03_sort"
MARKED_DATA_PATH="$GATK_PATH/01_marked"
RECALIBRATED_DATA_PATH="$GATK_PATH/02_BaseRecalibrator"
APPLYBQSR_DATA_PATH="$GATK_PATH/03_applyBQSR"
HAPLOTYPECALLER_DATA_PATH="$GATK_PATH/04_HaplotypeCaller"

# Set tools enviroment variables
PICARD_BIN="/mnt/raid6/bacphagenetwork/tools/picard_pre-built/v2.26.0/picard.jar"
GATK_BIN="/mnt/raid6/bacphagenetwork/tools/gatk-4.3.0.0/gatk"

# Prompt the paths
echo "The path to the genome data has been set to $GENOME_PATH."
echo "The indexing file is located at $INDEXING_FILE."
echo "The analysis path is set to $ANALYSIS_PATH."
echo "The samtools path is set to $SAMTOOLS_PATH."
echo "The GATK path is set to $GATK_PATH."
echo "The Picard binary is located at $PICARD_BIN."
echo "The GATK binary is located at $GATK_BIN."
echo "=============================="

# Check if the GATK output path exists
echo "Checking GATK output paths..."
if [ -d "$GATK_PATH" ]; then
    echo "GATK output path exists."
else
    echo "GATK output path does not exist. Creating directory..."
    mkdir -p "$GATK_PATH/01_marked"
    mkdir -p "$GATK_PATH/02_BaseRecalibrator"
    mkdir -p "$GATK_PATH/03_applyBQSR"
    mkdir -p "$GATK_PATH/04_HaplotypeCaller"
    echo "GATK output directory setup is complete." 
fi

# Check if the PICARD and GATK tools are available
echo "Checking if GATK and Picard tools are available..."
if [ ! -f "$PICARD_BIN" ]; then
    echo "Error: Picard tool not found at $PICARD_BIN."
    exit 1
fi
if [ ! -f "$GATK_BIN" ]; then
    echo "Error: GATK tool not found at $GATK_BIN."
    exit 1
fi
echo "GATK and Picard tools are available."
echo "=============================="

# Phase command line arguments
perform_mark_dulpicate=true
perform_index=true
perform_base_recalibrator=true
perform_apply_bqsr=true
perform_haplotype_caller=true

while getopts "mibah" opt; do
    case $opt in
        m) perform_mark_dulpicate=true ;;
        i) perform_index=true ;;
        b) perform_base_recalibrator=true ;;
        a) perform_apply_bqsr=true ;;
        h) perform_haplotype_caller=true ;;
        *) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
    esac
done

# If no options were provided, set all to true
if [ $OPTIND -eq 1 ]; then
    perform_mark_dulpicate=true
    perform_index=true
    perform_base_recalibrator=true
    perform_apply_bqsr=true
    perform_haplotype_caller=true
fi

# Get the list of all *.bam files
infile=($( cat sbatch.list | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print $0}' ))
sample_name=$(basename "$infile" .sorted.bam)
echo "Processing ${sample_name}..."

# Check if the sorted BAM file exists
if [ -f "$SORTED_DATA_PATH/${sample_name}.sorted.bam" ]; then
    echo "Sorted BAM file for ${sample_name} exists."
else
    echo "Error: Sorted BAM file for ${sample_name} does not exist."
    exit 1
fi

# Step 1: Use Picard to mark duplicates
if [ "$perform_mark_dulpicate" = true ]; then
    echo "Marking duplicates for ${sample_name}..."
    java -jar $PICARD_BIN MarkDuplicates \
        INPUT="$SORTED_DATA_PATH/${sample_name}.sorted.bam" \
        OUTPUT="$MARKED_DATA_PATH/${sample_name}.marked.bam" \
        METRICS_FILE="$MARKED_DATA_PATH/${sample_name}.marked.metrics.txt" \
        ASSUME_SORTED=true \
        REMOVE_DUPLICATES=false \
        VALIDATION_STRINGENCY=SILENT \
    || { echo "Error: Marking duplicates for $SORTED_DATA_PATH/${sample_name}.sorted.bam failed."; exit 1; }
    echo "Marking duplicates for ${sample_name} completed."
    echo "================================="
else
    echo "Skipping MarkDuplicates step for ${sample_name}."
fi

# Step 2: Index the marked BAM file
if [ "$perform_index" = true ]; then
    echo "Indexing the marked BAM file for ${sample_name}..."
    samtools index -@ 4 -b "$MARKED_DATA_PATH/${sample_name}.marked.bam" "$MARKED_DATA_PATH/${sample_name}.marked.bai" \
    || { echo "Error: Indexing the marked BAM file for ${sample_name} failed."; exit 1; }
    echo "Indexing the marked BAM file for ${sample_name} completed."
    echo "================================="
else
    echo "Skipping indexing step for ${sample_name}."
fi

# Step 3: Perform BaseRecalibrator
if [ "$perform_base_recalibrator" = true ]; then
    echo "Performing BaseRecalibrator for ${sample_name}..."
    $GATK_BIN BaseRecalibrator \
        -I "$MARKED_DATA_PATH/${sample_name}.marked.bam" \
        -R "$INDEXING_FILE" \
        --known-sites "$KNOWN_SITES_FILE" \
        -O "$RECALIBRATED_DATA_PATH/${sample_name}.recal_data.table" \
        --use-original-qualities \
    || { echo "Error: BaseRecalibrator for ${sample_name} failed"; exit 1; }
    echo "BaseRecalibrator for ${sample_name} completed."
    echo "================================="
else
    echo "Skipping BaseRecalibrator step for ${sample_name}."
fi

# Step 4: Apply BQSR
if [ "$perform_apply_bqsr" = true ]; then
    echo "Applying BQSR for ${sample_name}..."
    $GATK_BIN ApplyBQSR \
        -I "$MARKED_DATA_PATH/${sample_name}.marked.bam" \
        -R "$INDEXING_FILE" \
        --bqsr-recal-file "$RECALIBRATED_DATA_PATH/${sample_name}.recal_data.table" \
        -O "$APPLYBQSR_DATA_PATH/${sample_name}.recalibrated.bam" \
    || { echo "Error: ApplyBQSR for ${sample_name} failed"; exit 1; }
    echo "ApplyBQSR for ${sample_name} completed."
    echo "================================="
else
    echo "Skipping ApplyBQSR step for ${sample_name}."
fi

# Step 5: Perform HaplotypeCaller
if [ "$perform_haplotype_caller" = true ]; then
    echo "Performing HaplotypeCaller for ${sample_name}..."
    $GATK_BIN HaplotypeCaller \
        -I "$APPLYBQSR_DATA_PATH/${sample_name}.recalibrated.bam" \
        -R "$INDEXING_FILE" \
        -O "$HAPLOTYPECALLER_DATA_PATH/${sample_name}.g.vcf.gz" \
        -ERC GVCF \
        --native-pair-hmm-threads 2 \
    || { echo "Error: HaplotypeCaller for ${sample_name} failed"; exit 1; }
    echo "HaplotypeCaller for ${sample_name} completed."
    echo "================================="
else
    echo "Skipping HaplotypeCaller step for ${sample_name}."
fi

echo "All steps for ${sample_name} completed successfully."
echo "The output files are saved in the following directories:"
echo "1. $MARKED_DATA_PATH/${sample_name}.marked.bam"
echo "2. $MARKED_DATA_PATH/${sample_name}.marked.bai"
echo "3. $RECALIBRATED_DATA_PATH/${sample_name}.recal_data.table"
echo "4. $APPLYBQSR_DATA_PATH/${sample_name}.recalibrated.bam"
echo "5. $APPLYBQSR_DATA_PATH/${sample_name}.recalibrated.bai"
echo "6. $HAPLOTYPECALLER_DATA_PATH/${sample_name}.g.vcf.gz"
echo "7. $HAPLOTYPECALLER_DATA_PATH/${sample_name}.g.vcf.gz.tbi"
echo "================================="
echo "The analysis of ${sample_name} has been completed."
echo "All done."
