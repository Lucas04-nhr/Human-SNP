#!/bin/bash
#SBATCH --job-name=bwa_analysis
#SBATCH --output=./log/bwa_index.%j.out
#SBATCH --error=./log/bwa_index.%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=GENOME_PATH='/mnt/raid6/bacphagenetwork/data/skin_metagenome/Beijing/02_rm_host',INDEXING_PATH='/mnt/raid6/bacphagenetwork/data/bwa_index/chm13v2.0_noY.fa',ANALYSIS_PATH='/mnt/raid6/bacphagenetwork/data/bwa_analysis',INDEXING_FILE='/mnt/raid6/bacphagenetwork/data/bwa_index/chm13v2.0_noY.fa'

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

cd $INDEXING_PATH

echo "Current working directory: $(pwd)"

# Indexing
bwa index -a bwtsw $INDEXING_FILE