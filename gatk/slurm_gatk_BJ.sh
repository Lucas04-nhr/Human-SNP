#!/bin/bash

#SBATCH --job-name=gatk4_snp_calling
#SBATCH --output=./log/02/Beijing/gatk4_snp_calling_%j.out
#SBATCH --error=./log/02/Beijing/gatk4_snp_calling_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --export=ANALYSIS_PATH='/mnt/raid6/bacphagenetwork/data/samtools_analysis/Beijing',RESULTS_PATH='/mnt/raid6/bacphagenetwork/data/gatk4_analysis/Beijing'
#SBATCH --array=1-201%4

# Check whether the environment exists
if conda env list | grep -q "wescall"
then
    echo "Great! The environment already exists."
    # Activate the environment
    echo "Activating the environment..."
    conda activate wescall
else
    echo "Creating the environment..."
    conda env create -n wescall -f ../requirements.txt
    echo "The environment has been created, activating it..."
    conda activate wescall
fi

echo "Initialization is complete."

echo "The working directory has been changed to $ANALYSIS_PATH."

infile=($( cat bj_01_sbatch.list | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print $0}' ))
sample_name=$(echo "$infile" | grep -oE 'BJ[0-9]{3}')

# Check if the BAM file exists
if [ ! -f "$ANALYSIS_PATH/${sample_name}.bam" ]
then
    echo "Error: $ANALYSIS_PATH/${sample_name}.bam does not exist."
    exit 1
fi

# Run GATK4 HaplotypeCaller
echo "Processing $sample_name..."
echo "Running GATK4 HaplotypeCaller..."
echo "The path to the BAM file is $ANALYSIS_PATH/${sample_name}.bam."
echo "The path to the VCF file is $RESULTS_PATH/${sample_name}.vcf.gz."
gatk HaplotypeCaller \
   -R /path/to/reference.fasta \
   -I $ANALYSIS_PATH/${sample_name}.bam \
   -O $RESULTS_PATH/${sample_name}.vcf.gz \
   --native-pair-hmm-threads 4 || { echo "Error: GATK4 HaplotypeCaller failed in processing $sample_name."; exit 1; }

echo "The $sample_name BAM file has been successfully processed by GATK4 HaplotypeCaller."