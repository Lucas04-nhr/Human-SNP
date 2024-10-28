#!/bin/bash
#SBATCH --job-name=processing_BJ
#SBATCH --output=./log/Beijing/processing_BJ.%j.out
#SBATCH --error=./log/Beijing/processing_BJ.%j.err
#SBATCH --cpus-per-task=5
#SBATCH --mem=2G
#SBATCH --export=BASE_PATH='/mnt/raid6/bacphagenetwork/data/',GATK_OLD_BIN="/mnt/raid6/bacphagenetwork/tools/gatk-4.3.0.0/gatk",GATK_NEW_BIN="/mnt/raid6/bacphagenetwork/tools/gatk-4.5.0.0/gatk",JAVA_HOME='/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/',JAVA_BIN='/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/bin/java',LDFLAGS='-L/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/lib/server',CPPFLAGS='-I/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/include'
#SBATCH --array=1-201%4

# Initialize the environment
echo "Initializing..."


# Define input and output files
INPUT_VCF="$BASE_PATH/test/08_GenotypeGVCF/joint_genotyped.vcf.gz"
OUTPUT_VCF="$BASE_PATH/test/08_GenotypeGVCF/filtered_output.vcf"
REFERENCE_GENOME="$BASE_PATH/00_bwa_index/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa"

# Define filter expressions
SNP_FILTER="QD < 2.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"


# Apply hard filters to SNPs
gatk VariantFiltration \
    -R $REFERENCE_GENOME \
    -V $INPUT_VCF \
    --filter-expression "$SNP_FILTER" \
    --filter-name "SNP_Hard_Filter" \
    -O $OUTPUT_VCF

echo "Variant filtration complete. Filtered VCF is saved as $OUTPUT_VCF"