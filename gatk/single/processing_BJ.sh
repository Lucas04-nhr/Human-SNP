#!/bin/bash
#SBATCH --job-name=processing_BJ
#SBATCH --output=./log/Beijing/processing_BJ.%j.out
#SBATCH --error=./log/Beijing/processing_BJ.%j.err
#SBATCH --cpus-per-task=5
#SBATCH --mem=2G
#SBATCH --export=BASE_PATH='/mnt/raid6/bacphagenetwork/data/',GATK_OLD_BIN="/mnt/raid6/bacphagenetwork/tools/gatk-4.3.0.0/gatk",GATK_NEW_BIN="/mnt/raid6/bacphagenetwork/tools/gatk-4.5.0.0/gatk",JAVA_HOME='/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/',JAVA_BIN='/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/bin/java',LDFLAGS='-L/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/lib/server',CPPFLAGS='-I/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/include'
#SBATCH --array=1-8%8

# Initialize the environment
echo "Initializing..."

# Set the paths of the output files
export INDEXING_PATH="$BASE_PATH/00_bwa_index/GRCh38"
export INDEXING_FILE="$INDEXING_PATH/Homo_sapiens.GRCh38.dna.toplevel.fa"
export KNOWN_SITES_PATH="$BASE_PATH/00_bwa_index/GRCh38/known-sites/dbsnp138"
export KNOWN_SITES_FILE="$KNOWN_SITES_PATH/hg38_v0_Homo_sapiens_assembly38.dbsnp138.modified.vcf"
export SORTED_DATA_PATH="$BASE_PATH/03_sort/Beijing"
export RECALIBRATED_DATA_PATH="$BASE_PATH/05_BaseRecalibrator/Beijing"
export APPLYBQSR_DATA_PATH="$BASE_PATH/06_ApplyBQSR/Beijing"
export HAPLOTYPECALLER_DATA_PATH="$BASE_PATH/07_HaplotypeCaller/Beijing"

echo "The sorted *.bam files are located in $SORTED_DATA_PATH."
echo "The indexing data is located in $INDEXING_PATH."
echo "The indexing genome data is $INDEXING_FILE."
echo "The dbSNP database is located in $KNOWN_SITES_PATH."
echo "The dbSNP database file is $KNOWN_SITES_FILE."
echo "The recalibrated *.bam files will be saved in $RECALIBRATED_DATA_PATH."
echo "The ApplyBQSR results will be saved in $APPLYBQSR_DATA_PATH."
echo "The HaplotypeCaller results will be saved in $HAPLOTYPECALLER_DATA_PATH."

echo "Initializing completed."
echo "=============================="

# Get the list of all *.bam files
infile=($( cat BJ_sbatch.list | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print $0}' ))
sample_name=$(basename "$infile" .bam)

# Parse command line arguments
perform_base_recalibrator=false
perform_apply_bqsr=false
perform_haplotype_caller=false

while getopts "bah" opt; do
  case $opt in
    b) perform_base_recalibrator=true ;;
    a) perform_apply_bqsr=true ;;
    h) perform_haplotype_caller=true ;;
    *) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
  esac
done

# If no options were provided, set all to true
if ! $perform_base_recalibrator && ! $perform_apply_bqsr && ! $perform_haplotype_caller; then
  perform_base_recalibrator=true
  perform_apply_bqsr=true
  perform_haplotype_caller=true
fi

# Perform the BaseRecalibrator
if $perform_base_recalibrator; then
  echo "Performing BaseRecalibrator for ${sample_name}..."
  $GATK_OLD_BIN BaseRecalibrator \
    -I $SORTED_DATA_PATH/${sample_name}.bam \
    -R $INDEXING_FILE \
    --known-sites $KNOWN_SITES_FILE \
    -O $RECALIBRATED_DATA_PATH/${sample_name}.recal_data.table \
    --use-original-qualities \
  || { echo "BaseRecalibrator for ${sample_name} failed"; exit 1; }

  echo "BaseRecalibrator for ${sample_name} completed."
  echo "=============================="
else
  echo "BaseRecalibrator for ${sample_name} skipped."
  echo "=============================="
fi

# Perform the ApplyBQSR
if $perform_apply_bqsr; then
  echo "Performing ApplyBQSR for ${sample_name}..."
  $GATK_OLD_BIN ApplyBQSR \
    -I $SORTED_DATA_PATH/${sample_name}.bam \
    -R $INDEXING_FILE \
    --bqsr-recal-file $RECALIBRATED_DATA_PATH/${sample_name}.recal_data.table \
    -O $APPLYBQSR_DATA_PATH/${sample_name}.recalibrated.bam \
  || { echo "ApplyBQSR for ${sample_name} failed"; exit 1; }

  echo "ApplyBQSR for ${sample_name} completed."
  echo "=============================="
else
  echo "ApplyBQSR for ${sample_name} skipped."
  echo "=============================="
fi

# Perform the HaplotypeCaller
if $perform_haplotype_caller; then
  echo "Performing HaplotypeCaller for ${sample_name}..."
  $GATK_OLD_BIN HaplotypeCaller \
    -I $APPLYBQSR_DATA_PATH/${sample_name}.recalibrated.bam \
    -R $INDEXING_FILE \
    -O $HAPLOTYPECALLER_DATA_PATH/${sample_name}.g.vcf.gz \
    -ERC GVCF \
    --native-pair-hmm-threads 2 \
  || { echo "HaplotypeCaller for ${sample_name} failed"; exit 1; }

  echo "HaplotypeCaller for ${sample_name} completed."
  echo "=============================="
else
  echo "HaplotypeCaller for ${sample_name} skipped."
  echo "=============================="
fi

echo "All processes completed with no error.(exit 0)"
