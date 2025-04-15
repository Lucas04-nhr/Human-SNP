#! /bin/bash
#SBATCH --job-name=ninglab_data_gatk_joint
#SBATCH --output=./ninglab_data_gatk_joint.%j.out
#SBATCH --error=./ninglab_data_gatk_joint.%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G

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
INDEXING_FILE='/mnt/raid6/bacphagenetwork/data/00_bwa_index/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa'
ANALYSIS_PATH="${BASE_PATH}/01_bwa_analysis"
SAMTOOLS_PATH="${BASE_PATH}/02_samtools"
GATK_PATH="${BASE_PATH}/03_gatk"

# Set known-sites path
KNOWN_SITES_BASE_PATH="$BASE_PATH/../00_bwa_index/GRCh38/known-sites"
KNOWN_SITES_1000G="$KNOWN_SITES_BASE_PATH/1000g/hg38_v0_1000G_phase1.snps.high_confidence.hg38.modified.vcf"
KNOWN_SITES_DBSNP="$KNOWN_SITES_BASE_PATH/dbsnp138/hg38_v0_Homo_sapiens_assembly38.dbsnp138.modified.vcf"
KNOWN_SITES_HAPMAP="$KNOWN_SITES_BASE_PATH/hapmap/hg38_v0_hapmap_3.3.hg38.modified.vcf"
KNOWN_SITES_OMNI="$KNOWN_SITES_BASE_PATH/omni/hg38_v0_1000G_omni2.5.hg38.modified.vcf"

# Set other environment variables
KNOWN_SITES_PATH="$BASE_PATH/../00_bwa_index/GRCh38/known-sites/dbsnp138"
KNOWN_SITES_FILE="$KNOWN_SITES_PATH/hg38_v0_Homo_sapiens_assembly38.dbsnp138.modified.vcf"
SORTED_DATA_PATH="$SAMTOOLS_PATH/03_sort"
MARKED_DATA_PATH="$GATK_PATH/01_marked"
RECALIBRATED_DATA_PATH="$GATK_PATH/02_BaseRecalibrator"
APPLYBQSR_DATA_PATH="$GATK_PATH/03_applyBQSR"
HAPLOTYPECALLER_DATA_PATH="$GATK_PATH/04_HaplotypeCaller"
GVCF_DATA_PATH="$GATK_PATH/05_GenotypeGVCF"
VARIANTRECALIBRATOR_DATA_PATH="$GATK_PATH/06_VariantRecalibrator"
APPLYVQSR_DATA_PATH="$GATK_PATH/07_ApplyVQSR"

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

# Check if the known-site files exist
echo "Checking known-site files..."
if [ -f "$KNOWN_SITES_1000G" ]; then
    echo "1000G known-site file exists."
else
    echo "1000G known-site file does not exist. Please check the path."
    exit 1
fi
if [ -f "$KNOWN_SITES_DBSNP" ]; then
    echo "dbSNP known-site file exists."
else
    echo "dbSNP known-site file does not exist. Please check the path."
    exit 1
fi
if [ -f "$KNOWN_SITES_HAPMAP" ]; then
    echo "HapMap known-site file exists."
else
    echo "HapMap known-site file does not exist. Please check the path."
    exit 1
fi
if [ -f "$KNOWN_SITES_OMNI" ]; then
    echo "Omni known-site file exists."
else
    echo "Omni known-site file does not exist. Please check the path."
    exit 1
fi
echo "All known-site files exist."
echo "=============================="

# Check if the GATK JOINT CALL output paths exist
echo "Checking GVCF output path..."
if [ -d "$GVCF_DATA_PATH" ]; then
    echo "GVCF output path exists."
else
    echo "GVCF output path does not exist. Creating directory..."
    mkdir -p "$GVCF_DATA_PATH"
    echo "GVCF output directory setup is complete." 
fi

echo "Checking VariantRecalibrator output path..."
if [ -d "$VARIANTRECALIBRATOR_DATA_PATH" ]; then
    echo "VariantRecalibrator output path exists."
else
    echo "VariantRecalibrator output path does not exist. Creating directory..."
    mkdir -p "$VARIANTRECALIBRATOR_DATA_PATH"
    echo "VariantRecalibrator output directory setup is complete." 
fi

echo "Checking ApplyVQSR output path..."
if [ -d "$APPLYVQSR_DATA_PATH" ]; then
    echo "ApplyVQSR output path exists."
else
    echo "ApplyVQSR output path does not exist. Creating directory..."
    mkdir -p "$APPLYVQSR_DATA_PATH"
    echo "ApplyVQSR output directory setup is complete." 
fi
echo "=============================="

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
perform_merge=false
perform_genotype_gvcf=false
perform_variant_recalibrator=false
perform_apply_vqsr=false

while getopts "mgvr" opt; do
  case $opt in
    m) perform_merge=true ;;
    g) perform_genotype_gvcf=true ;;
    v) perform_variant_recalibrator=true ;;
    r) perform_apply_vqsr=true ;;
    *) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
  esac
done

# If no options were provided, set all to true
if ! $perform_merge && ! $perform_genotype_gvcf && ! $perform_variant_recalibrator && ! $perform_apply_vqsr; then
  perform_merge=true
  perform_genotype_gvcf=true
  perform_variant_recalibrator=true
  perform_apply_vqsr=true
fi

# Merge all gVCF files
if $perform_merge; then
  # Get the list of all *.g.vcf files
  GVCF_FILES=("$HAPLOTYPECALLER_DATA_PATH"/*.g.vcf.gz)

  gvcf_list=$(printf " -V %s" "${GVCF_FILES[@]}")

  echo "Merging all gVCF files..."
  $GATK_BIN CombineGVCFs \
    -R $INDEXING_FILE \
    $gvcf_list \
    -O $GVCF_DATA_PATH/merged_genome.vcf.gz \
  || { echo "Merging gVCF files failed"; exit 1; }

  echo "Merging gVCF files completed."
  echo "=============================="
else
  echo "Skipping to merge all gVCFs..."
  echo "=============================="
fi

# Perform joint genotyping
if $perform_genotype_gvcf; then
  echo "Performing joint genotyping..."
  $GATK_BIN GenotypeGVCFs \
    -R $INDEXING_FILE \
    -V $GVCF_DATA_PATH/merged_genome.vcf.gz \
    -O $GVCF_DATA_PATH/joint_genotyped.vcf.gz \
  || { echo "Joint genotyping failed"; exit 1; }

  echo "Joint genotyping completed."
  echo "=============================="
else
  echo "Skipping joint genotyping..."
  echo "=============================="
fi

# Perform VariantRecalibrator
if $perform_variant_recalibrator; then
  echo "Performing VariantRecalibrator..."
  $GATK_BIN VariantRecalibrator \
    -R $INDEXING_FILE \
    -V $GVCF_DATA_PATH/joint_genotyped.vcf.gz \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $KNOWN_SITES_HAPMAP \
    --resource:omni,known=false,training=true,truth=true,prior=12.0 $KNOWN_SITES_OMNI \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 $KNOWN_SITES_1000G \
    --resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $KNOWN_SITES_DBSNP \
    -mode SNP \
    -an QD -an FS -an MQRankSum -an ReadPosRankSum -an SOR \
    -O $VARIANTRECALIBRATOR_DATA_PATH/joint_genotyped.recal \
    --tranches-file $VARIANTRECALIBRATOR_DATA_PATH/joint_genotyped.tranches \
    --rscript-file $VARIANTRECALIBRATOR_DATA_PATH/joint_genotyped.plots.R \
  || { echo "VariantRecalibrator failed"; exit 1; }

  echo "VariantRecalibrator completed."
  echo "=============================="
else
  echo "Skipping VariantRecalibrator..."
  echo "=============================="
fi

# Perform ApplyVQSR
if $perform_apply_vqsr; then
  echo "Performing ApplyVQSR..."
  $GATK_BIN ApplyVQSR \
    -R $INDEXING_FILE \
    -V $GVCF_DATA_PATH/joint_genotyped.vcf.gz \
    --recal-file $VARIANTRECALIBRATOR_DATA_PATH/joint_genotyped.recal \
    --tranches-file $VARIANTRECALIBRATOR_DATA_PATH/joint_genotyped.tranches \
    -mode SNP \
    -O $APPLYVQSR_DATA_PATH/joint_genotyped.vcf.gz \
  || { echo "ApplyVQSR failed"; exit 1; }

  echo "ApplyVQSR completed."
  echo "=============================="
else
  echo "Skipping ApplyVQSR..."
  echo "=============================="
fi

# Print the final message
echo "All steps completed successfully."
echo "The output files are saved in the following directories:"
echo "1. $GVCF_DATA_PATH/merged_genome.vcf.gz"
echo "2. $GVCF_DATA_PATH/joint_genotyped.vcf.gz"
echo "3. $VARIANTRECALIBRATOR_DATA_PATH/joint_genotyped.recal"
echo "4. $VARIANTRECALIBRATOR_DATA_PATH/joint_genotyped.tranches"
echo "5. $VARIANTRECALIBRATOR_DATA_PATH/joint_genotyped.plots.R"
echo "6. $APPLYVQSR_DATA_PATH/joint_genotyped.vcf.gz"
echo "7. $APPLYVQSR_DATA_PATH/joint_genotyped.vcf.gz.tbi"
echo "=============================="
echo "The analysis has been completed."
