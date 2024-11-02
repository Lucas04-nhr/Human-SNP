#! /bin/bash
#SBATCH --job-name=plink_BJ
#SBATCH --output=./BJ_log.%j.out
#SBATCH --error=./BJ_log.%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

# Initialize the environment
echo "Initializing the environment..."
echo "=============================="
echo ""


# Load conda
echo "Activating conda environment..."
source /home/bacphagenetwork/.bashrc
conda activate snp_analysis
echo "The conda environment has been activated."
echo "=============================="

# Set the paths of the output files
echo "Setting the paths of the output files..."
export JAVA_HOME='/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1'
export JAVA_BIN='/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/bin/java'
export LDFLAGS='-L/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/lib/server'
export CPPFLAGS='-I/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/include'
export BASE_PATH="/mnt/raid6/bacphagenetwork/data"

export GATK_OLD_BIN="/mnt/raid6/bacphagenetwork/tools/gatk-4.3.0.0/gatk"
export GATK_NEW_BIN="/mnt/raid6/bacphagenetwork/tools/gatk-4.5.0.0/gatk"
export PICARD_OLD_BIN='/mnt/raid6/bacphagenetwork/tools/picard_pre-built/v2.26.0/picard.jar'
export PICARD_NEW_BIN='/mnt/raid6/bacphagenetwork/tools/picard_pre-built/v3.0/picard.jar'
export PLINK_OLD_BIN="/mnt/raid6/bacphagenetwork/tools/plink-1.07-x86_64/plink"
export PLINK_NEW_BIN="/mnt/raid6/bacphagenetwork/tools/plink_1.9_linux_x86_64/plink"

export INDEXING_PATH="$BASE_PATH/00_bwa_index/GRCh38"
export INDEXING_FILE="$INDEXING_PATH/Homo_sapiens.GRCh38.dna.toplevel.fa"

export KNOWN_SITES_BASE_PATH="$BASE_PATH/00_bwa_index/GRCh38/known-sites"
export KNOWN_SITES_1000G="$KNOWN_SITES_BASE_PATH/1000g/hg38_v0_1000G_phase1.snps.high_confidence.hg38.modified.vcf"
export KNOWN_SITES_DBSNP="$KNOWN_SITES_BASE_PATH/dbsnp138/hg38_v0_Homo_sapiens_assembly38.dbsnp138.modified.vcf"
export KNOWN_SITES_HAPMAP="$KNOWN_SITES_BASE_PATH/hapmap/hg38_v0_hapmap_3.3.hg38.modified.vcf"
export KNOWN_SITES_OMNI="$KNOWN_SITES_BASE_PATH/omni/hg38_v0_1000G_omni2.5.hg38.modified.vcf"

export UNFILTERED_GVCF_PATH="$BASE_PATH/08_GenotypeGVCF/Beijing"
export FILTERED_GVCF_PATH="$BASE_PATH/10_ApplyVQSR/Beijing"
export PLINK_PATH="$BASE_PATH/12_plink/Beijing"

echo "The UNFILTERED GenotypeGVCF results is located in $UNFILTERED_GVCF_PATH."
echo "The FILTERED GenotypeGVCF result is located in $FILTERED_GVCFPATH."
echo "The plink files will be located in $PLINK_PATH."
# Add prompts of more sub-folders of plink here...

echo "Setting completed."

echo ""

echo "Initializing completed."
echo "=============================="

# Performing plink converting
echo "Converting the VCF files to plink format..."
$PLINK_NEW_BIN --noweb --vcf $FILETERED_GVCF_PATH/joint_genotyped.filtered.vcf.gz --recode --allow-extra-chr --out $PLINK_PATH/converted_genotyped \
|| { echo "Error: plink converting failed."; exit 1; }

echo "The plink converting has been completed."
echo "=============================="

# Performing plink 
echo "Performing..."
$PLINK_NEW_BIN --noweb --file $PLINK_PATH/converted_genotyped --out $PLINK_PATH/converted_genotyped --allow-extra-chr --make-bed \
|| { echo "Error: plink executing failed."; exit 1; }

echo "The process has been completed."
