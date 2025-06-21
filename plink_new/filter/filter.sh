#! /bin/bash
#SBATCH --job-name=plink_volcano_plot
#SBATCH --output=./volcano_plot_log.%j.out
#SBATCH --error=./volcano_plot_log.%j.err
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
export PLINK_CONVERT_PATH="$PLINK_PATH/converted"
export PLINK_FILTER_RESULT="$PLINK_PATH/filtered"

echo "The UNFILTERED GenotypeGVCF results is located in $UNFILTERED_GVCF_PATH."
echo "The FILTERED GenotypeGVCF result is located in $FILTERED_GVCF_PATH."
echo "The plink result files with bacteria added are located in $PLINK_CONVERT_PATH."
echo "The volcano plot files will be located in $PLINK_FILTER_RESULT."

mkdir -p $PLINK_FILTER_RESULT

echo "Setting completed."

echo ""

echo "Initializing completed."
echo "=============================="

# Execute the script

for file in $PLINK_CONVERT_PATH/*.assoc.linear; do
  echo "Processing file: $(basename $file)"
  python plink_result_filter.py --input-file $file --output-directory $PLINK_FILTER_RESULT || {
    echo "An error occurred while running the script on $(basename $file)."
    exit 1
  }
  echo "=============================="
  echo ""
done
