#!/bin/bash
#SBATCH --job-name=convert_test
#SBATCH --output=./log/convert_gvcf_test.%j.out
#SBATCH --error=./log/convert_gvcf_test.%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --export=BASE_PATH='/mnt/raid6/bacphagenetwork/data/',GATK_OLD_BIN="/mnt/raid6/bacphagenetwork/tools/gatk-4.3.0.0/gatk",GATK_NEW_BIN="/mnt/raid6/bacphagenetwork/tools/gatk-4.5.0.0/gatk",JAVA_HOME='/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/',JAVA_BIN='/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/bin/java',LDFLAGS='-L/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/lib/server',CPPFLAGS='-I/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/include'
#SBATCH --array=1-88%3

# Initialize the environment
echo "Initializing..."

# Activate conda env
source /home/bacphagenetwork/.bashrc
conda activate snp_analysis

# Set the paths of the output files
export JAVA_HOME='/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/'
export JAVA_BIN='/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/bin/java'
export LDFLAGS='-L/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/lib/server'
export CPPFLAGS='-I/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/include'
export BASE_PATH="/mnt/raid6/bacphagenetwork/data"

export INDEXING_PATH="$BASE_PATH/00_bwa_index/GRCh38"
export INDEXING_FILE="$INDEXING_PATH/Homo_sapiens.GRCh38.dna.toplevel.fa"

export GATK_OLD_BIN="/mnt/raid6/bacphagenetwork/tools/gatk-4.3.0.0/gatk"
export GATK_NEW_BIN="/mnt/raid6/bacphagenetwork/tools/gatk-4.5.0.0/gatk"
export BCFTOOLS_BIN="/mnt/raid6/bacphagenetwork/tools/bcftools/bin/bcftools"

export HAPLOTYPECALLER_DATA_PATH="$BASE_PATH/07_HaplotypeCaller/Beijing"
export GVCF_TO_VCF_PATH="$BASE_PATH/11_ConvertedVCF/Beijing"
export TEST_VCF_PATH="$BASE_PATH/test/11_ConvertedVCF"

echo "The index file is located in $INDEXING_PATH."
echo "The HaplotypeCaller results is located in $HAPLOTYPECALLER_DATA_PATH."
echo "The GVCF to VCF results is located in $GVCF_TO_VCF_PATH."
echo "The test VCF results will be simlinked to $TEST_VCF_PATH."

echo "Initializing completed."
echo "=============================="

# Get the list of all *.g.vcf.gz files
infile=($( cat test_sbatch.list | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print $0}' ))
sample_name=$(basename "$infile" .g.vcf.gz)

# Perform Converting
echo "Converting ${sample_name}..."
$BCFTOOLS_BIN convert --gvcf2vcf $infile --fasta-ref $INDEXING_FILE -O z -o $GVCF_TO_VCF_PATH/$sample_name.vcf.gz \
|| { echo "Converting ${sample_name} failed."; exit 1; }

# Perform Simlinking
echo "Simlinking ${sample_name}..."
ln -s $GVCF_TO_VCF_PATH/$sample_name.vcf.gz $TEST_VCF_PATH/$sample_name.vcf.gz \
|| { echo "Simlinking ${sample_name} failed."; exit 1; }

echo "Converting ${sample_name} completed."
