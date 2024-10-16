#! /bin/bash
# Initialize the environment
echo "Initializing..."

# Set the paths of the output files
export JAVA_HOME='/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/'
export JAVA_BIN='/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/bin/java'
export LDFLAGS='-L/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/lib/server'
export CPPFLAGS='-I/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/include'
export BASE_PATH="/mnt/raid6/bacphagenetwork/data/test"

export GATK_OLD_BIN="/mnt/raid6/bacphagenetwork/tools/gatk-4.3.0.0/gatk"
export GATK_NEW_BIN="/mnt/raid6/bacphagenetwork/tools/gatk-4.5.0.0/gatk"

export INDEXING_PATH="$BASE_PATH/../00_bwa_index/GRCh38"
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

echo "The sorted *.bam files are located in $SORTED_DATA_PATH."
echo "The indexing data is located in $INDEXING_PATH."
echo "The indexing genome data is $INDEXING_FILE."
echo "The recalibrated *.bam files is located in $RECALIBRATED_DATA_PATH."
echo "The ApplyBQSR results is located in $APPLYBQSR_DATA_PATH."
echo "The HaplotypeCaller results is located in $HAPLOTYPECALLER_DATA_PATH."
echo "The GenotypeGVCF results will be located in $GENOTYPE_GVCF_PATH."
echo "The VariantRecalibrator results will be located in $VARIANTRECALIBRATOR_DATA_PATH."
echo "The ApplyVQSR results will be located in $APPLYVQSR_DATA_PATH."

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
  gvcf_files=($(ls $HAPLOTYPECALLER_DATA_PATH/*.g.vcf.gz))
  gvcf_list=$(printf " -V %s" "${gvcf_files[@]}")

  echo "Merging all gVCF files..."
  $GATK_OLD_BIN GenomicsDBImport \
    -R $INDEXING_FILE \
    --genomicsdb-workspace-path $GENOTYPE_GVCF_PATH/genomicsdb \
    --batch-size 20 \
    $gvcf_list \
    --overwrite-existing-genomicsdb-workspace true \
    -L 1 -L 2 -L 3 -L 4 -L 5 -L 6 -L 7 -L 8 -L 9 \
    -L 10 -L 11 -L 12 -L 13 -L 14 -L 15 -L 16 -L 17 \
    -L 18 -L 19 -L 20 -L 21 -L 22 -L X \
    --max-num-intervals-to-import-in-parallel 2 \
    --reader-threads 5 \
  || { echo "GenomicsDBImport failed"; exit 1; }

  $GATK_OLD_BIN GenotypeGVCFs \
    -R $INDEXING_FILE \
    -V gendb://$GENOTYPE_GVCF_PATH/genomicsdb \
    -O $GENOTYPE_GVCF_PATH/merged_genome.vcf.gz \
    --overwrite-existing-genomicsdb-workspace true \
    -L 1 -L 2 -L 3 -L 4 -L 5 -L 6 -L 7 -L 8 -L 9 \
    -L 10 -L 11 -L 12 -L 13 -L 14 -L 15 -L 16 -L 17 \
    -L 18 -L 19 -L 20 -L 21 -L 22 -L X \
    --max-num-intervals-to-import-in-parallel 2 \
    --reader-threads 5 \
  || { echo "Genotyping from GenomicsDB failed"; exit 1; }

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
    -V $GENOTYPE_GVCF_PATH/merged_genome.vcf.gz \
    -O $GENOTYPE_GVCF_PATH/joint_genotyped.vcf.gz \
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
    -V $GENOTYPE_GVCF_PATH/joint_genotyped.vcf.gz \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $KNOWN_SITES_HAPMAP \
    --resource:omni,known=false,training=true,truth=true,prior=12.0 $KNOWN_SITES_OMNI \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 $KNOWN_SITES_1000G \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $KNOWN_SITES_DBSNP \
    -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum \
    -mode SNP \
    -O $VARIANTRECALIBRATOR_DATA_PATH/output.recal \
    --tranches-file $VARIANTRECALIBRATOR_DATA_PATH/output.tranches \
    --rscript-file $VARIANTRECALIBRATOR_DATA_PATH/output.plots.R \
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
    -V $GENOTYPE_GVCF_PATH/joint_genotyped.vcf.gz \
    -O $APPLYVQSR_DATA_PATH/joint_genotyped.filtered.vcf.gz \
    --recal-file $VARIANTRECALIBRATOR_DATA_PATH/output.recal \
    --tranches-file $VARIANTRECALIBRATOR_DATA_PATH/output.tranches \
    -mode SNP \
  || { echo "Applying VQSR to joint genotyped VCF failed"; exit 1; }

  echo "VQSR completed."
  echo "=============================="
else
  echo "Skipping VQSR for joint genotyped VCF."
  echo "=============================="
fi

echo "All processes completed."
