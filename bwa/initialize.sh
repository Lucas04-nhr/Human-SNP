#! /bin/bash

# Check whether conda is installed
if ! command -v conda &> /dev/null
then
    echo "Conda could not be found, please install it first."
    exit
fi

# Check whether the environment exists
if conda env list | grep -q "bwa"
then
    echo "Great! The environment already exists."
    # Activate the environment
    echo "Activating the environment..."
    conda activate bwa
else
    echo "Creating the environment..."
    conda env create -f requirements.txt
    echo "The environment has been created."
fi

echo "Initialization is complete."

# Set the path to the genome data
echo "Please enter the path to the genome data:"
read genome_path
# The path of the genome data is '/mnt/raid6/bacphagenetwork/data/skin_metagenome/Beijing/02_rm_host'
export GENOME_PATH=$genome_path
echo "The path to the genome data has been set to $GENOME_PATH."

# Set the path to the indexing data
echo "Please enter where the the indexing data located:"
read indexing_path

export INDEXING_PATH=$indexing_path
echo "The path to the indexing data has been set to $INDEXING_PATH."


# Indexing
bwa index -a bwtsw $INDEXING_PATH


# Analyse the genome data
echo "Analysing the genome data..."
for file_fq1 in $(ls *_1.fastq.gz)
do
    file_name=$(basename "$file_fq1")
    sample_name="${file_name%%_1*}"
    file_fq2="${sample_name}_2.fastq.gz"
    bwa index -p $sample_name $GENOME_PATH/$sample_name.fasta
done