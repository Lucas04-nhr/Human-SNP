#! /bin/bash

export PATH=/mnt/raid6/bacphagenetwork/tools/bcftools/bin/:$PATH

BASE_PATH="/mnt/raid6/bacphagenetwork/data/00_bwa_index/GRCh38/known-sites/"

# Define arrays of input and output file paths
declare -A VCF_FILES=(
    ["$BASE_PATH/1000g/.original/hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf"]="$BASE_PATH/1000g/hg38_v0_1000G_phase1.snps.high_confidence.hg38.modified.vcf"
    ["$BASE_PATH/hapmap/.original/hg38_v0_hapmap_3.3.hg38.vcf"]="$BASE_PATH/hapmap/hg38_v0_hapmap_3.3.hg38.modified.vcf"
    ["$BASE_PATH/omni/.original/hg38_v0_1000G_omni2.5.hg38.vcf"]="$BASE_PATH/omni/hg38_v0_1000G_omni2.5.hg38.modified.vcf"
)

# Loop through each set of input and output paths
for INPUT_VCF in "${!VCF_FILES[@]}"; do
    OUTPUT_VCF="${VCF_FILES[$INPUT_VCF]}"
    OUTPUT_INDEX="${OUTPUT_VCF}.idx"

    echo "Modifying chromosome names in $INPUT_VCF..."

    # Modify chromosome names in the VCF file
    bcftools annotate --rename-chrs ./chr_names.txt "$INPUT_VCF" -o "$OUTPUT_VCF"

    echo "=============================="
    echo "Reindexing $OUTPUT_VCF..."

    # Reindex new VCF file
    gatk IndexFeatureFile -I "$OUTPUT_VCF"

    echo "Chromosome names have been modified and saved to $OUTPUT_VCF."
    echo "Output file has already been indexed and saved to $OUTPUT_INDEX."
    echo "=============================="
done

echo "All chromosome names have been modified and reindexed successfully."
