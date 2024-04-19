#! /bin/bash

# Set the path to the genome data
echo "Please enter the path to the genome data:"
read genome_path
# The path of the genome data is '/mnt/raid6/bacphagenetwork/data/skin_metagenome/Beijing/02_rm_host'
export GENOME_PATH=$genome_path
echo "The path to the genome data has been set to $GENOME_PATH."

# Check whether the index exists
if [ ! -f "./samples.index" ]
then
    echo "The index does not exist, creating it ..."

else
    echo "The index already exists, do you want to overwrite it? (y/n)"
    read answer
    if [ $answer == "y" ]
    then
        echo "Overwriting the index ..."
    else
        echo "The index will not be overwritten."
        exit
    fi
fi

# Add the genome data to the config file "samples.index"
echo "Adding the genome data to the config file..."
for file_fq1 in $(ls ${GENOME_PATH}/*_1.fastq.gz)
do
    file_name=$(basename "$file_fq1")
    sample_name="${file_name%%_1*}"
    echo "Processing $sample_name..."
    file_fq2="${sample_name}_2.fastq.gz"
    echo "$sample_name  $file_fq1   0" >> samples.index
    echo "The genome data of $sample_name has been added to the config file."
done