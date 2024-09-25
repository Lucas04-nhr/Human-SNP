#!/bin/bash
#SBATCH --job-name=processing_BJ
#SBATCH --output=./log/Beijing/processing_BJ.%j.out
#SBATCH --error=./log/Beijing/processing_BJ.%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=BASE_PATH='/mnt/raid6/bacphagenetwork/data/',INDEXING_PATH='/mnt/raid6/bacphagenetwork/data/00_bwa_index/chm13v2/chm13v2.0_noY.fa',GATK_OLD_BIN="/mnt/raid6/bacphagenetwork/tools/gatk-4.3.0.0/gatk",GATK_NEW_BIN="/mnt/raid6/bacphagenetwork/tools/gatk-4.5.0.0/gatk",
#SBATCH --array=1-201%4

# Initialize the environment
echo "Initializing..."

# Set the paths of the output files
INDEXING_FILE="$BASE_PATH/00_bwa_index/chm13v2/chm13v2.0_noY.fa"
KNOWN_SITES_PATH="$BASE_PATH/00_bwa_index/chm13v2/known_sites"
SORTED_DATA_PATH="$BASE_PATH/03_sort/Beijing"
RECALIBRATED_DATA_PATH="$BASE_PATH/05_BaseRecalibrator/Beijing"
APPLYBQSR_DATA_PATH="$BASE_PATH/06_ApplyBQSR/Beijing"
HAPLOTYPECALLER_DATA_PATH="$BASE_PATH/07_HaplotypeCaller/Beijing"

echo "The sorted *.bam files are located in $SORTED_DATA_PATH."
echo "The recalibrated *.bam files will be saved in $RECALIBRATED_DATA_PATH."
echo "The ApplyBQSR results will be saved in $APPLYBQSR_DATA_PATH."
echo "The HaplotypeCaller results will be saved in $HAPLOTYPECALLER_DATA_PATH."

echo "Initializing completed."
echo "=============================="

# Get the list of all *.bam files
infile=($( cat BJ_sbatch.list | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print $0}' ))
sample_name=$(basename "$infile" .sorted.bam)

# Perform the BaseRecalibrator
echo "Performing BaseRecalibrator for ${sample_name}..."
$GATK_NEW_BIN BaseRecalibrator \
  -I $SORTED_DATA_PATH/${sample_name}.sorted.bam \
  -R $INDEXING_PATH \
  --known-sites $KNOWN_SITES_PATH/dbsnp_146.hg38.vcf \

