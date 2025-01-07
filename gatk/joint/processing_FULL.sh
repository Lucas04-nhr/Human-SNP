#! /bin/bash
#SBATCH --job-name=joint_full
#SBATCH --output=./log_full.%j.out
#SBATCH --error=./log_full.%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --export=BASE_PATH='/mnt/raid6/bacphagenetwork/data/',GATK_OLD_BIN="/mnt/raid6/bacphagenetwork/tools/gatk-4.3.0.0/gatk",GATK_NEW_BIN="/mnt/raid6/bacphagenetwork/tools/gatk-4.5.0.0/gatk",JAVA_HOME='/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/',JAVA_BIN='/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/bin/java',LDFLAGS='-L/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/lib/server',CPPFLAGS='-I/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/include'

# Initialize the environment
echo "Initializing..."

# Activate conda env
source /home/bacphagenetwork/.bashrc
conda activate snp_analysis

# Set the paths of the output files
export JAVA_HOME='/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1'
export JAVA_BIN='/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/bin/java'
export LDFLAGS='-L/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/lib/server'
export CPPFLAGS='-I/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/include'
export BASE_PATH="/mnt/raid6/bacphagenetwork/data"

export GATK_OLD_BIN="/mnt/raid6/bacphagenetwork/tools/gatk-4.3.0.0/gatk"
export GATK_NEW_BIN="/mnt/raid6/bacphagenetwork/tools/gatk-4.5.0.0/gatk"

export INDEXING_PATH="$BASE_PATH/00_bwa_index/GRCh38"
export INDEXING_FILE="$INDEXING_PATH/Homo_sapiens.GRCh38.dna.toplevel.fa"

export KNOWN_SITES_BASE_PATH="$BASE_PATH/00_bwa_index/GRCh38/known-sites"
export KNOWN_SITES_1000G="$KNOWN_SITES_BASE_PATH/1000g/hg38_v0_1000G_phase1.snps.high_confidence.hg38.modified.vcf"
export KNOWN_SITES_DBSNP="$KNOWN_SITES_BASE_PATH/dbsnp138/hg38_v0_Homo_sapiens_assembly38.dbsnp138.modified.vcf"
export KNOWN_SITES_HAPMAP="$KNOWN_SITES_BASE_PATH/hapmap/hg38_v0_hapmap_3.3.hg38.modified.vcf"
export KNOWN_SITES_OMNI="$KNOWN_SITES_BASE_PATH/omni/hg38_v0_1000G_omni2.5.hg38.modified.vcf"

export SORTED_DATA_PATH_BJ="$BASE_PATH/03_sort/Beijing"
export RECALIBRATED_DATA_PATH_BJ="$BASE_PATH/05_BaseRecalibrator/Beijing"
export APPLYBQSR_DATA_PATH_BJ="$BASE_PATH/06_ApplyBQSR/Beijing"
export HAPLOTYPECALLER_DATA_PATH_BJ="$BASE_PATH/07_HaplotypeCaller/Beijing"
export GENOTYPE_GVCF_PATH_BJ="$BASE_PATH/08_GenotypeGVCF/Beijing"
export VARIANTRECALIBRATOR_DATA_PATH_BJ="$BASE_PATH/09_VariantRecalibrator/Beijing"
export APPLYVQSR_DATA_PATH_BJ="$BASE_PATH/10_ApplyVQSR/Beijing"

export SORTED_DATA_PATH_GZ="$BASE_PATH/03_sort/Guangzhou"
export RECALIBRATED_DATA_PATH_GZ="$BASE_PATH/05_BaseRecalibrator/Guangzhou"
export APPLYBQSR_DATA_PATH_GZ="$BASE_PATH/06_ApplyBQSR/Guangzhou"
export HAPLOTYPECALLER_DATA_PATH_GZ="$BASE_PATH/07_HaplotypeCaller/Guangzhou"
export GENOTYPE_GVCF_PATH_GZ="$BASE_PATH/08_GenotypeGVCF/Guangzhou"
export VARIANTRECALIBRATOR_DATA_PATH_GZ="$BASE_PATH/09_VariantRecalibrator/Guangzhou"
export APPLYVQSR_DATA_PATH_GZ="$BASE_PATH/10_ApplyVQSR/Guangzhou"

export GENOTYPE_GVCF_PATH_FULL="$BASE_PATH/08_GenotypeGVCF/Full"
export VARIANTRECALIBRATOR_DATA_PATH_FULL="$BASE_PATH/09_VariantRecalibrator/Full"
export APPLYVQSR_DATA_PATH_FULL="$BASE_PATH/10_ApplyVQSR/Full"



echo "The sorted *.bam files are located in $SORTED_DATA_PATH_BJ and $SORTED_DATA_PATH_GZ."
echo "The indexing data is located in $INDEXING_PATH_BJ and $INDEXING_PATH_GZ."
echo "The indexing genome data is $INDEXING_FILE_BJ and $INDEXING_FILE_GZ."
echo "The recalibrated *.bam files is located in $RECALIBRATED_DATA_PATH_BJ and $RECALIBRATED_DATA_PATH_GZ."
echo "The ApplyBQSR results is located in $APPLYBQSR_DATA_PATH_BJ and $APPLYBQSR_DATA_PATH_GZ."
echo "The HaplotypeCaller results is located in $HAPLOTYPECALLER_DATA_PATH_BJ and $HAPLOTYPECALLER_DATA_PATH_GZ."
echo "The GenotypeGVCF results will be located in $GENOTYPE_GVCF_PATH_FULL."
echo "The VariantRecalibrator results will be located in $VARIANTRECALIBRATOR_DATA_PATH_FULL."
echo "The ApplyVQSR results will be located in $APPLYVQSR_DATA_PATH_FULL."

echo "Initializing completed."
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

if $perform_merge; then
  # Merge all gVCF files
  gvcf_files_bj=($(ls $HAPLOTYPECALLER_DATA_PATH_BJ/*.g.vcf.gz))
  gvcf_files_gz=($(ls $HAPLOTYPECALLER_DATA_PATH_GZ/*.g.vcf.gz))
  gvcf_files_full=("${gvcf_files_bj[@]}" "${gvcf_files_gz[@]}")

  gvcf_list=$(printf " -V %s" "${gvcf_files_full[@]}")

  echo "Merging all gVCF files..."
  $GATK_OLD_BIN CombineGVCFs \
    -R $INDEXING_FILE \
    $gvcf_list \
    -O $GENOTYPE_GVCF_PATH_FULL/merged_genome.vcf.gz \
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
  $GATK_OLD_BIN GenotypeGVCFs \
    -R $INDEXING_FILE \
    -V $GENOTYPE_GVCF_PATH_FULL/merged_genome.vcf.gz \
    -O $GENOTYPE_GVCF_PATH_FULL/joint_genotyped.vcf.gz \
  || { echo "Joint genotyping failed"; exit 1; }

  echo "Joint genotyping completed."
  echo "=============================="
else
  echo "Skipping joint genotyping..."
  echo "=============================="
fi

# Perform VQSR
if $perform_variant_recalibrator; then
  echo "Performing VQSR for joint genotyped VCF..."
  $GATK_OLD_BIN VariantRecalibrator \
    -R $INDEXING_FILE \
    -V $GENOTYPE_GVCF_PATH_FULL/joint_genotyped.vcf.gz \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $KNOWN_SITES_HAPMAP \
    --resource:omni,known=false,training=true,truth=true,prior=12.0 $KNOWN_SITES_OMNI \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 $KNOWN_SITES_1000G \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $KNOWN_SITES_DBSNP \
    -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum \
    -mode SNP \
    -O $VARIANTRECALIBRATOR_DATA_PATH_FULL/output.recal \
    --tranches-file $VARIANTRECALIBRATOR_DATA_PATH_FULL/output.tranches \
    --rscript-file $VARIANTRECALIBRATOR_DATA_PATH_FULL/output.plots.R \
  || { echo "VQSR for joint genotyped VCF failed"; exit 1; }

  echo "VQSR for joint genotyped VCF completed."
  echo "=============================="
else
  echo "Skipping VQSR for joint genotyped VCF."
  echo "=============================="
fi

if $perform_apply_vqsr; then
  echo "Applying VQSR to joint genotyped VCF..."
  $GATK_OLD_BIN ApplyVQSR \
    -R $INDEXING_FILE \
    -V $GENOTYPE_GVCF_PATH_FULL/joint_genotyped.vcf.gz \
    -O $APPLYVQSR_DATA_PATH_FULL/joint_genotyped.filtered.vcf.gz \
    --recal-file $VARIANTRECALIBRATOR_DATA_PATH_FULL/output.recal \
    --tranches-file $VARIANTRECALIBRATOR_DATA_PATH_FULL/output.tranches \
    -mode SNP \
  || { echo "Applying VQSR to joint genotyped VCF failed"; exit 1; }

  echo "VQSR completed."
  echo "=============================="
else
  echo "Skipping VQSR for joint genotyped VCF."
  echo "=============================="
fi

echo "All processes completed."
