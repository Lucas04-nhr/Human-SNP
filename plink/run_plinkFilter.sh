#!/bin/bash
#SBATCH --job-name=plink_filter_job  # 作业名称
#SBATCH --output=plink_filter.out  # 标准输出和错误日志文件
#SBATCH --error=plink_filter.err
#SBATCH --mem=128G  # 内存需求

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

export UNFILTERED_GVCF_PATH="$BASE_PATH/08_GenotypeGVCF/Full"
export FILTERED_GVCF_PATH="$BASE_PATH/10_ApplyVQSR/Full"
export PLINK_PATH="$BASE_PATH/12_plink/Full"
export PLINK_OUTPUT_PATH="$PLINK_PATH/output/bac_age"
export PLINK_CORRECTION_PATH="$PLINK_PATH/correction"
export PLINK_RESULT_PATH="$PLINK_PATH/results/bac_age"
export PLINK_MODIFIED_PATH="$PLINK_PATH/modified"
export PLINK_MERGE_PATH="$PLINK_PATH/merged"
  
echo "The UNFILTERED GenotypeGVCF results is located in $UNFILTERED_GVCF_PATH."
echo "The FILTERED GenotypeGVCF result is located in $FILTERED_GVCF_PATH."
echo "The plink output files will be located in $PLINK_OUTPUT_PATH."
echo "The plink correction files will be located in $PLINK_CORRECTION_PATH."
echo "The plink result files will be located in $PLINK_RESULT_PATH."
echo "The plink modified files will be located in $PLINK_MODIFIED_PATH."
# Add prompts of more sub-folders of plink here...

# Run the plinkFilter.py
python plinkFilter.py --directory $PLINK_MODIFIED_PATH/bac_age --output_file $PLINK_MERGE_PATH/filtered_bac_age.csv || \
{ echo "The plinkFilter.py has failed." ; exit 1; }
echo "The plinkFilter.py has been executed successfully."