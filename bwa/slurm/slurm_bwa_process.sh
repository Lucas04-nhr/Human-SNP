#!/bin/bash
#SBATCH --job-name=bwa_analysis
#SBATCH --output=bwa_analysis.out
#SBATCH --error=bwa_analysis.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=GENOME_PATH='/mnt/raid6/bacphagenetwork/data/skin_metagenome/Beijing/02_rm_host',INDEXING_PATH='/mnt/raid6/bacphagenetwork/data/bwa_index/chm13v2.0_noY.fa',ANALYSIS_PATH='/mnt/raid6/bacphagenetwork/data/bwa_analysis'

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

# # The path of the genome data is '/mnt/raid6/bacphagenetwork/data/skin_metagenome/Beijing/02_rm_host'
echo "The path to the genome data has been set to $GENOME_PATH."

# # The path to the indexing data is '/mnt/raid6/bacphagenetwork/data/bwa_index/chm13v2.0_noY.fa'
echo "The path to the indexing data has been set to $INDEXING_PATH."

# # The path to store the analysis results is '/mnt/raid6/bacphagenetwork/data/bwa_analysis'
echo "The path to store the analysis results has been set to $ANALYSIS_PATH."

# Indexing
bwa index -a bwtsw $INDEXING_PATH

