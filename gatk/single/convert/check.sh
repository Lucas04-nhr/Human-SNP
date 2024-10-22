#! /bin/bash

# Initialize the environment
echo "Initializing..."

# Set the paths of the output files
export JAVA_HOME='/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/'
export JAVA_BIN='/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/bin/java'
export LDFLAGS='-L/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/lib/server'
export CPPFLAGS='-I/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/include'
export BASE_PATH="/mnt/raid6/bacphagenetwork/data/"

export GATK_OLD_BIN="/mnt/raid6/bacphagenetwork/tools/gatk-4.3.0.0/gatk"
export GATK_NEW_BIN="/mnt/raid6/bacphagenetwork/tools/gatk-4.5.0.0/gatk"

export INDEXING_PATH="$BASE_PATH/00_bwa_index/GRCh38"
export INDEXING_FILE="$INDEXING_PATH/Homo_sapiens.GRCh38.dna.toplevel.fa"

export KNOWN_SITES_BASE_PATH="$BASE_PATH/00_bwa_index/GRCh38/known-sites"
export KNOWN_SITES_1000G="$KNOWN_SITES_BASE_PATH/1000g/hg38_v0_1000G_phase1.snps.high_confidence.modified.hg38.vcf"
export KNOWN_SITES_DBSNP="$KNOWN_SITES_BASE_PATH/dbsnp138/hg38_v0_Homo_sapiens_assembly38.dbsnp138.modified.vcf"
export KNOWN_SITES_HAPMAP="$KNOWN_SITES_BASE_PATH/hapmap/hg38_v0_hapmap_3.3.hg38.modified.vcf"
export KNOWN_SITES_OMNI="$KNOWN_SITES_BASE_PATH/omni/hg38_v0_1000G_omni2.5.hg38.modified.vcf"

export SORTED_DATA_PATH="$BASE_PATH/03_sort"
export RECALIBRATED_DATA_PATH="$BASE_PATH/05_BaseRecalibrator"
export APPLYBQSR_DATA_PATH="$BASE_PATH/06_ApplyBQSR"
export HAPLOTYPECALLER_DATA_PATH="$BASE_PATH/07_HaplotypeCaller"
export GENOTYPE_GVCF_PATH="$BASE_PATH/08_GenotypeGVCF"
export VARIANTRECALIBRATOR_DATA_PATH="$BASE_PATH/09_VariantRecalibrator"
export APPLYVQSR_DATA_PATH="$BASE_PATH/10_ApplyVQSR"
export GVCF_TO_VCF_PATH="$BASE_PATH/11_ConvertedVCF"

echo "The sorted *.bam files are located in $SORTED_DATA_PATH."
echo "The indexing data is located in $INDEXING_PATH."
echo "The indexing genome data is $INDEXING_FILE."
echo "The recalibrated *.bam files is located in $RECALIBRATED_DATA_PATH."
echo "The ApplyBQSR results is located in $APPLYBQSR_DATA_PATH."
echo "The HaplotypeCaller results is located in $HAPLOTYPECALLER_DATA_PATH."
echo "The GenotypeGVCF results will be located in $GENOTYPE_GVCF_PATH."
echo "The VariantRecalibrator results will be located in $VARIANTRECALIBRATOR_DATA_PATH."
echo "The ApplyVQSR results will be located in $APPLYVQSR_DATA_PATH."

# Check if the input directory exists

# Function to check if a directory exists
check_directory() {
  if [ ! -d "$1" ]; then
    echo "Directory $1 does not exist. Please check manually."
    exit 1
  else
    echo "Directory $1 exists."
  fi
}

# Function to check if a file exists
check_file() {
  if [ ! -f "$1" ]; then
    echo "File $1 does not exist. Please check manually."
    exit 1
  else
    echo "File $1 exists."
  fi
}

# Check directories
check_directory "$SORTED_DATA_PATH"
check_directory "$INDEXING_PATH"
check_directory "$KNOWN_SITES_BASE_PATH"
check_directory "$RECALIBRATED_DATA_PATH"
check_directory "$APPLYBQSR_DATA_PATH"
check_directory "$HAPLOTYPECALLER_DATA_PATH"

# Check files
check_file "$INDEXING_FILE"
check_file "$KNOWN_SITES_DBSNP"
check_file "$KNOWN_SITES_HAPMAP"
check_file "$KNOWN_SITES_OMNI"


# Check if the output directory exists

# Function to check and create directories if they do not exist
check_and_create_directory() {
  if [ ! -d "$1" ]; then
    echo "Directory $1 does not exist. Creating it..."
    mkdir -p "$1/Beijing" "$1/Guangzhou"
  else
    echo "Directory $1 exists."
  fi
}

# Check and create necessary directories
check_and_create_directory "$GENOTYPE_GVCF_PATH"
check_and_create_directory "$VARIANTRECALIBRATOR_DATA_PATH"
check_and_create_directory "$APPLYVQSR_DATA_PATH"

echo "Initializing completed."
