#! /bin/bash

export PATH=/mnt/raid6/bacphagenetwork/tools/bcftools/bin/:$PATH

# Define the input and output file paths
INPUT_VCF="/mnt/raid6/bacphagenetwork/data/00_bwa_index/GRCh38/known-sites/Homo_sapiens_assembly38.dbsnp138.vcf"
OUTPUT_VCF="/mnt/raid6/bacphagenetwork/data/00_bwa_index/GRCh38/known-sites/Homo_sapiens_assembly38_modified.dbsnp138.vcf"
OUTPUT_INDEX="/mnt/raid6/bacphagenetwork/data/00_bwa_index/GRCh38/known-sites/Homo_sapiens_assembly38_modified.dbsnp138.idx"

# Modify chromosome names in the VCF file
bcftools annotate --rename-chrs ./chr_names.txt $INPUT_VCF -o $OUTPUT_VCF

# Reindex new VCF file
bcftools index $OUTPUT_VCF -f $OUTPUT_INDEX

echo "Chromosome names have been modified and saved to $OUTPUT_VCF."
echo "Output file has already been indexed and saved to $OUTPUT_INDEX."
