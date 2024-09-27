#!/bin/bash

# Initialize the environment
echo "Initializing..."

# Set the paths of the output files
BASE_PATH="/mnt/raid6/bacphagenetwork/data/"
GATK_OLD_BIN="/mnt/raid6/bacphagenetwork/tools/gatk-4.3.0.0/gatk"
GATK_NEW_BIN="/mnt/raid6/bacphagenetwork/tools/gatk-4.5.0.0/gatk"
INDEXING_PATH="$OUTPUT_BASE_PATH/00_bwa_index/GRCh38"
INDEXING_FILE="$INDEXING_PATH/Homo_sapiens.GRCh38.dna.toplevel.fa"
KNOWN_SITES_PATH="$BASE_PATH/00_bwa_index/GRCh38/known_sites"
KNOWN_SITES_FILE="$KNOWN_SITES_PATH/Homo_sapiens_assembly38.dbsnp138.vcf"
SORTED_DATA_PATH="$BASE_PATH/03_sort/Beijing"
RECALIBRATED_DATA_PATH="$BASE_PATH/05_BaseRecalibrator/Beijing"
APPLYBQSR_DATA_PATH="$BASE_PATH/06_ApplyBQSR/Beijing"
HAPLOTYPECALLER_DATA_PATH="$BASE_PATH/07_HaplotypeCaller/Beijing"

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
infile="/mnt/raid6/bacphagenetwork/data/test/BJ001.bam"
sample_name=$(basename "$infile" .bam)

# Perform the BaseRecalibrator
echo "Performing BaseRecalibrator for ${sample_name}..."
$GATK_NEW_BIN BaseRecalibrator \
  -I $SORTED_DATA_PATH/${sample_name}.bam \
  -R $INDEXING_FILE \
  --known-sites $KNOWN_SITES_FILE \
  -O $RECALIBRATED_DATA_PATH/${sample_name}.recal_data.table \
  --use-original-qualities \
  --num-threads 5 \
|| { echo "BaseRecalibrator for ${sample_name} failed"; exit 1; }

echo "BaseRecalibrator for ${sample_name} completed."
echo "=============================="

# Perform the ApplyBQSR

echo "Performing ApplyBQSR for ${sample_name}..."
$GATK_NEW_BIN ApplyBQSR \
  -I $SORTED_DATA_PATH/${sample_name}.bam \
  -R $INDEXING_PATH \
  --bqsr-recal-file $RECALIBRATED_DATA_PATH/${sample_name}.recal_data.table \
  -O $APPLYBQSR_DATA_PATH/${sample_name}.recalibrated.bam \
  --num-threads 5 \
|| { echo "ApplyBQSR for ${sample_name} failed"; exit 1; }

echo "ApplyBQSR for ${sample_name} completed."
echo "=============================="

# Perform the HaplotypeCaller

echo "Performing HaplotypeCaller for ${sample_name}..."
$GATK_NEW_BIN HaplotypeCaller \
  -I $APPLYBQSR_DATA_PATH/${sample_name}.recalibrated.bam \
  -R $INDEXING_PATH \
  -O $HAPLOTYPECALLER_DATA_PATH/${sample_name}.g.vcf.gz \
  --native-pair-hmm-threads 5 \
|| { echo "HaplotypeCaller for ${sample_name} failed"; exit 1; }

echo "HaplotypeCaller for ${sample_name} completed."
echo "=============================="

echo "All processes completed."
