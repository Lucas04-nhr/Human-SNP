#! /bin/bash
#SBATCH --job-name=plink_full
#SBATCH --output=./full_log.%j.out
#SBATCH --error=./full_log.%j.err
#SBATCH --cpus-per-task=20
#SBATCH --mem=256G

# Initialize the environment
echo "Initializing the environment..."
echo "=============================="
echo ""


# Load conda
echo "Activating conda environment..."
source /home/bacphagenetwork/.bashrc
source /mnt/raid3/bacphagenetwork/miniconda3/bin/activate
conda init
conda activate snp_analysis
echo "The conda environment has been activated."
echo "=============================="

# Set the paths of the output files
echo "Setting the paths of the output files..."
export BASE_PATH="/mnt/raid6/bacphagenetwork/data"

export INDEXING_PATH="$BASE_PATH/00_bwa_index/GRCh38"
export INDEXING_FILE="$INDEXING_PATH/Homo_sapiens.GRCh38.dna.toplevel.fa"

export KNOWN_SITES_BASE_PATH="$BASE_PATH/00_bwa_index/GRCh38/known-sites"
export KNOWN_SITES_1000G="$KNOWN_SITES_BASE_PATH/1000g/hg38_v0_1000G_phase1.snps.high_confidence.hg38.modified.vcf"
export KNOWN_SITES_DBSNP="$KNOWN_SITES_BASE_PATH/dbsnp138/hg38_v0_Homo_sapiens_assembly38.dbsnp138.modified.vcf"
export KNOWN_SITES_HAPMAP="$KNOWN_SITES_BASE_PATH/hapmap/hg38_v0_hapmap_3.3.hg38.modified.vcf"
export KNOWN_SITES_OMNI="$KNOWN_SITES_BASE_PATH/omni/hg38_v0_1000G_omni2.5.hg38.modified.vcf"

export UNFILTERED_GVCF_PATH="$BASE_PATH/08_GenotypeGVCF/Full"
export FILTERED_GVCF_PATH="$BASE_PATH/10_ApplyVQSR/Full"
export PLINK_PATH="$BASE_PATH/12_plink_Full"
export PLINK_OUTPUT_PATH="$PLINK_PATH/output"
export PLINK_RESULT_PATH="$PLINK_PATH/results"

echo "The UNFILTERED GenotypeGVCF results is located in $UNFILTERED_GVCF_PATH."
echo "The FILTERED GenotypeGVCF result is located in $FILTERED_GVCF_PATH."
echo "The plink output files will be located in $PLINK_OUTPUT_PATH."
echo "The plink result files will be located in $PLINK_RESULT_PATH."
# Add prompts of more sub-folders of plink here...

mkdir -p $PLINK_OUTPUT_PATH
mkdir -p $PLINK_RESULT_PATH

echo "Setting completed."

echo ""

echo "Initializing completed."
echo "=============================="

# Phase command arguments

plink_convert=false
plink_preprocess=false
plink_execute=false
plink_draw_fig=false
covar_number=""

# sbatch plink_full.sh -e --covar-number=4

while getopts "cpe-:" opt; do
  case $opt in
    c) plink_convert=true ;;
    p) plink_preprocess=true ;;
    e) plink_execute=true ;;
    -)
      case "${OPTARG}" in
        covar-number=*)
          covar_number=${OPTARG#*=}
          ;;
        *)
          echo "Invalid option: --${OPTARG}" ;;
      esac
      ;;
    \?) echo "Invalid option: $OPTARG" ;;
  esac
done

# If no options were provided, set all to true

if ! $plink_convert && ! $plink_preprocess && ! $plink_execute; then
  plink_convert=true
  plink_preprocess=true
  plink_execute=true
fi

# Performing plink converting
if $plink_convert; then
  echo "Converting the VCF files to plink format..."
  plink --noweb --vcf $FILTERED_GVCF_PATH/joint_genotyped.filtered.vcf.gz --recode --allow-extra-chr --out $PLINK_OUTPUT_PATH/converted_genotyped \
  || { echo "Error: plink converting failed."; exit 1; }
  echo "The plink converting has been completed."
  echo "=============================="
else
  echo "The plink converting has been skipped."
  echo "=============================="
fi

# Performing plink preprocessing
if $plink_preprocess; then
  echo "Performing plink preprocessing..."
  plink --noweb --file $PLINK_OUTPUT_PATH/converted_genotyped --set-missing-var-ids @:# --recode --out $PLINK_OUTPUT_PATH/converted_genotyped --allow-extra-chr --make-bed \
  || { echo "Error: plink preprocessing failed."; exit 1; }
  echo "The plink preprocessing has been completed."
  echo "=============================="
else
  echo "The plink preprocessing has been skipped."
  echo "=============================="
fi

# Performing plink execution
if $plink_execute; then
  if [ "$covar_number" -gt 483 ] || [ "$covar_number" -lt 3 ]; then
    echo "Error: covar_number must be a number between 3 and 483, set to 5 by default."
    covar_number=5
  elif [ -z "$covar_number" ]; then
    echo "No covar_number provided, set to 5 by default."
    covar_number=5
  fi
  echo "The covar_number is set to $covar_number."

  echo "=============================="
  echo "Performing plink execution..."
  plink --bfile $PLINK_OUTPUT_PATH/converted_genotyped --linear --assoc --pheno $PLINK_PATH/phenotype_full.tsv --all-pheno --covar $PLINK_PATH/covariate_full.tsv --covar-number $covar_number --missing --out $PLINK_RESULT_PATH/result --noweb --allow-extra-chr --allow-no-sex \
  || { echo "Error: plink execution failed."; exit 1; }
  echo "The plink execution has been completed."
  echo "=============================="

else
  echo "The plink execution has been skipped."
  echo "=============================="
fi

echo "The process has been completed."
