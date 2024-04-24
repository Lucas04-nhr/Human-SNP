#! /bin/bash

# Check whether conda is installed
if ! command -v conda &> /dev/null
then
    echo "Conda could not be found, please install it first."
    exit
fi

conda init bash
source ~/.bashrc
conda activate base
# Check whether the environment exists
if conda env list | grep -q "bwa"
then
    echo "Great! The environment already exists."
    # Activate the environment
    echo "Activating the environment..."
    conda activate bwa
else
    echo "Creating the environment..."
    conda env create -n bwa -f requirements.txt
    echo "The environment has been created, activating it..."
    conda activate bwa
fi
export INDEXING_PATH='/mnt/raid6/bacphagenetwork/data/bwa_index'
export INDEXING_FILE='chm13v2.0_noY.fa'

cd $INDEXING_PATH

echo "Current working directory: $(pwd)"

# Indexing
bwa index -a bwtsw $INDEXING_FILE