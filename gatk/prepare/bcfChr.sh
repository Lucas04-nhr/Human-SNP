#!/bin/bash
#SBATCH --job-name=modify_vcf
#SBATCH --output=modify_vcf_%j.log
#SBATCH --error=modify_vcf_%j.err
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --mem=4G

# Load necessary modules
module load bcftools

# Define the input and output file paths
INPUT_VCF="/mnt/raid6/bacphagenetwork/data/00_bwa_index/GRCh38/known-sites/Homo_sapiens_assembly38.dbsnp138.vcf"
OUTPUT_VCF="/mnt/raid6/bacphagenetwork/data/00_bwa_index/GRCh38/known-sites/Homo_sapiens_assembly38_modified.dbsnp138.vcf"

# Modify chromosome names in the VCF file
bcftools annotate --rename-chrs chr_names.txt $INPUT_VCF -o $OUTPUT_VCF

echo "Chromosome names have been modified and saved to $OUTPUT_VCF"