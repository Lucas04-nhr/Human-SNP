export BCFTOOLS_BIN=“/mnt/raid6/bacphagenetwork/tools/bcftools/bin/bcftools”

# Define the input and output file paths
export INPUT_VCF="/mnt/raid6/bacphagenetwork/data/00_bwa_index/GRCh38/known-sites/Homo_sapiens_assembly38.dbsnp138.vcf"
export OUTPUT_VCF="/mnt/raid6/bacphagenetwork/data/00_bwa_index/GRCh38/known-sites/Homo_sapiens_assembly38_modified.dbsnp138.vcf"
export OUTPUT_INDEX="/mnt/raid6/bacphagenetwork/data/00_bwa_index/GRCh38/known-sites/Homo_sapiens_assembly38_modified.dbsnp138.idx"

# Modify chromosome names in the VCF file
$BCFTOOLS_BIN annotate --rename-chrs ./chr_names.txt $INPUT_VCF -o $OUTPUT_VCF

# Reindex new VCF file
$BCFTOOLS_BIN index $OUTPUT_VCF -f $OUTPUT_INDEX

echo "Chromosome names have been modified and saved to $OUTPUT_VCF."
echo "Output file has already been indexed and saved to $OUTPUT_INDEX."
