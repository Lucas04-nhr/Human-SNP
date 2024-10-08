#!/bin/bash
#SBATCH --job-name=bwa_analysis
#SBATCH --output=./log/Guangzhou/bwa_analysis.%j.out
#SBATCH --error=./log/Guangzhou/bwa_analysis.%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=GENOME_PATH='/mnt/raid6/bacphagenetwork/data/skin_metagenome/Guangzhou/02_rm_host',INDEXING_PATH='/mnt/raid6/bacphagenetwork/data/bwa_index/Homo_sapiens.GRCh38.dna.toplevel.fa',ANALYSIS_PATH='/mnt/raid6/bacphagenetwork/data/bwa_analysis/Guangzhou'
#SBATCH --array=1-160%4

conda init bash
source ~/.bashrc

# Check whether the environment exists
if conda env list | grep -q "bwa"
then
    echo "Great! The environment already exists."
    # Activate the environment
    echo "Activating the environment..."
    conda activate bwa
else
    echo "Creating the environment..."
    conda env create -n bwa -f ../requirements.txt
    echo "The environment has been created, activating it..."
    conda activate bwa
fi

echo "Initialization is complete."

# The path of the genome data is '/mnt/raid6/bacphagenetwork/data/skin_metagenome/Beijing/02_rm_host'
echo "The path to the genome data has been set to $GENOME_PATH."

# The path to the indexing data is '/mnt/raid6/bacphagenetwork/data/bwa_index/chm13v2.0_noY.fa'
echo "The path to the indexing data has been set to $INDEXING_PATH."

# The path to store the analysis results is '/mnt/raid6/bacphagenetwork/data/bwa_analysis'
echo "The path to store the analysis results has been set to $ANALYSIS_PATH."

# Indexing
# bwa index -a bwtsw $INDEXING_PATH

# Analyse the genome data
echo "Analysing the genome data..."
infile=($( cat gz_sbatch.list | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print $0}' ))
sample_name_temp=$(echo "$infile" | awk -F'/' '{n=split($NF,a,"_"); printf "%03d", a[n-1]}')
sample_name="GZ$sample_name_temp"

## TO DO
## Fix the file_fq1 and file_fq2 paths
## DONE

file_name=$(echo "$infile" | sed 's/_1.fastq.gz$//')

echo "Processing $sample_name..."
file_fq1="${file_name}_1.fastq.gz"
file_fq2="${file_name}_2.fastq.gz"
touch $ANALYSIS_PATH/${sample_name}.sam
bwa mem -t 4 $INDEXING_PATH $file_fq1 $file_fq2 -a > $ANALYSIS_PATH/log/${sample_name}.log > $ANALYSIS_PATH/${sample_name}.sam || { echo "Error: bwa mem failed in processing $sample_name."; exit 1; }
echo "The analysis of $sample_name has been completed."
