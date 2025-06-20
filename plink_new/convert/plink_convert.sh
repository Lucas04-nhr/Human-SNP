#! /bin/bash
#SBATCH --job-name=plink_convert
#SBATCH --output=./convert_log.%j.out
#SBATCH --error=./convert_log.%j.err
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
export PLINK_RESULT_PATH="$PLINK_PATH/results"
export PLINK_CONVERT_PATH="$PLINK_PATH/converted"

echo "The UNFILTERED GenotypeGVCF results is located in $UNFILTERED_GVCF_PATH."
echo "The FILTERED GenotypeGVCF result is located in $FILTERED_GVCF_PATH."
echo "The plink result files is located in $PLINK_RESULT_PATH."
echo "The plink result files with bacteria added will be located in $PLINK_CONVERT_PATH."
# Add prompts of more sub-folders of plink here...

mkdir -p $PLINK_CONVERT_PATH

echo "Setting completed."

echo ""

echo "Initializing completed."
echo "=============================="

# Execute the script

python ./plinkConvert.py --input-directory $PLINK_RESULT_PATH --output-directory $PLINK_CONVERT_PATH --pheno-file $PLINK_PATH/phenotype_full.tsv || {
    echo "An error occurred while running the script."
    exit 1
}
