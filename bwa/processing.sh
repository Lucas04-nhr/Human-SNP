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
    conda env create -n bwa -f requirements.txt
    echo "The environment has been created, activating it..."
    conda activate bwa
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

# Set the path to store analysis results
echo "Please enter the path to store the analysis results:"
read analysis_path
# The path to store the analysis results is '/mnt/raid6/bacphagenetwork/data/bwa_analysis'
export ANALYSIS_PATH=$analysis_path
echo "The path to store the analysis results has been set to $ANALYSIS_PATH."

# Indexing
bwa index -a bwtsw $INDEXING_PATH


# Analyse the genome data
echo "Analysing the genome data..."
for file_fq1 in $(ls *_1.fastq.gz)
do
    file_name=$(basename "$file_fq1")
    sample_name="${file_name%%_1*}"
    echo "Processing $sample_name..."
    file_fq2="${sample_name}_2.fastq.gz"
    bwa mem -t 4 -R "\@RG\tID:{library}\tLB:{library}\tPL:Illumina\tPU:{sample}\tSM:{sample}\" -p $INDEXING_PATH $file_fq1 $file_fq2 -a > $ANALYSIS_PATH/${sample_name}.sam > $ANALYSIS_PATH/log/${sample_name}.log || { echo "Error: bwa mem failed in processing $sample_name."; exit 1; }
    echo "The analysis of $sample_name has been completed."
done