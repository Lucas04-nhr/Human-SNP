#! /bin/bash

# Initialize the environment
echo "Initializing..."

# Set the paths of the output files
export JAVA_HOME='/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1'
export JAVA_BIN='/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/bin/java'
export LDFLAGS='-L/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/lib/server'
export CPPFLAGS='-I/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/include'
export BASE_PATH="/mnt/raid6/bacphagenetwork/data"
export TEST_PATH="$BASE_PATH/test"

export GATK_OLD_BIN="/mnt/raid6/bacphagenetwork/tools/gatk-4.3.0.0/gatk"
export GATK_NEW_BIN="/mnt/raid6/bacphagenetwork/tools/gatk-4.5.0.0/gatk"
export PLINK_BIN="/mnt/raid6/bacphagenetwork/tools/plink-1.07-x86_64/plink"

export INDEXING_PATH="$BASE_PATH/00_bwa_index/GRCh38"
export INDEXING_FILE="$INDEXING_PATH/Homo_sapiens.GRCh38.dna.toplevel.fa"

export KNOWN_SITES_BASE_PATH="$BASE_PATH/../00_bwa_index/GRCh38/known-sites"
export KNOWN_SITES_1000G="$KNOWN_SITES_BASE_PATH/1000g/hg38_v0_1000G_phase1.snps.high_confidence.hg38.modified.vcf"
export KNOWN_SITES_DBSNP="$KNOWN_SITES_BASE_PATH/dbsnp138/hg38_v0_Homo_sapiens_assembly38.dbsnp138.modified.vcf"
export KNOWN_SITES_HAPMAP="$KNOWN_SITES_BASE_PATH/hapmap/hg38_v0_hapmap_3.3.hg38.modified.vcf"
export KNOWN_SITES_OMNI="$KNOWN_SITES_BASE_PATH/omni/hg38_v0_1000G_omni2.5.hg38.modified.vcf"

export SORTED_DATA_PATH="$BASE_PATH"
export RECALIBRATED_DATA_PATH="$BASE_PATH"
export APPLYBQSR_DATA_PATH="$BASE_PATH"
export HAPLOTYPECALLER_DATA_PATH="$BASE_PATH/07_HaplotypeCaller"
export GENOTYPE_GVCF_PATH="$BASE_PATH/08_GenotypeGVCF"
export VARIANTRECALIBRATOR_DATA_PATH="$BASE_PATH/09_VariantRecalibrator"
export APPLYVQSR_DATA_PATH="$BASE_PATH/10_ApplyVQSR"
export TEST_VCF_PATH="$BASE_PATH/11_ConvertedVCF"
export PLINK_BASE_PATH="$BASE_PATH/12_plink"
export PLINK_TEST_PATH="$TEST_PATH/12_plink"
# Add more sub-folders of plink here...
export PLINK_CONVERTED_DATA="$PLINK_BASE_PATH/01_Converted"
export PLINK_CONVERTED_DATA_TEST="$PLINK_TEST_PATH/01_Converted"

echo "The sorted *.bam files are located in $SORTED_DATA_PATH."
echo "The indexing data is located in $INDEXING_PATH."
echo "The indexing genome data is $INDEXING_FILE."
echo "The recalibrated *.bam files are located in $RECALIBRATED_DATA_PATH."
echo "The ApplyBQSR results are located in $APPLYBQSR_DATA_PATH."
echo "The HaplotypeCaller results are located in $HAPLOTYPECALLER_DATA_PATH."
echo "The GenotypeGVCF results are located in $GENOTYPE_GVCF_PATH."
echo "The VariantRecalibrator results are located in $VARIANTRECALIBRATOR_DATA_PATH."
echo "The ApplyVQSR results are located in $APPLYVQSR_DATA_PATH."
echo "The converted VCF files are located in $TEST_VCF_PATH."
echo "The plink files will be located in $PLINK_BASE_PATH."
echo "The plink files for testing will be simlinked to $PLINK_TEST_PATH."
# Add prompts of more sub-folders of plink here...
echo "The converted plink files will be located in $PLINK_CONVERTED_DATA."
echo "The converted plink files for testing will be simlinked to $PLINK_CONVERTED_DATA_TEST."

# Function to check if a directory exists
check_directory() {
  if [ ! -d "$1" ]; then
    echo "Directory $1 does not exist. Please check manually."
    exit 1
  else
    echo "Directory $1 exists."
  fi
}

# Function to check if an output directory exists
check_output_directory() {
  if [ ! -d "$1" ]; then
    echo "Directory $1 does not exist. Creating it..."
    mkdir -p "$1"
  else
    echo "Directory $1 exists, would you like to overwrite it? (y/n)"
    read answer
    if [ $answer == "y" ] || [ $answer == "Y" ]; then
      echo "Overwriting the directory..."
      rm -rf "$1"
      mkdir -p "$1"
    else
      echo "The directory will not be overwritten."
    fi
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
check_directory "$GVCF_TO_VCF_PATH"
check_directory "$VARIANTRECALIBRATOR_DATA_PATH"
check_directory "$APPLYVQSR_DATA_PATH"
check_directory "$TEST_VCF_PATH"

# Check files
check_file "$JAVA_BIN"
check_file "$GATK_OLD_BIN"
check_file "$GATK_NEW_BIN"
check_file "$PLINK_BIN"

check_file "$INDEXING_FILE"
check_file "$KNOWN_SITES_DBSNP"
check_file "$KNOWN_SITES_HAPMAP"
check_file "$KNOWN_SITES_OMNI"
check_file "$KNOWN_SITES_1000G"

# Check if the output directory exists
check_output_directory "$PLINK_BASE_PATH"
check_output_directory "$PLINK_TEST_PATH"
check_output_directory "$PLINK_CONVERTED_DATA"
check_output_directory "$PLINK_CONVERTED_DATA_TEST"


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
