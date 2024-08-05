#! /bin/bash

# Set the path to the genome data
echo "Please enter the path to the pre-processed *.bam data of Beijing:"
read genome_path
# The path of the pre-processed *.bam data is '/mnt/raid6/bacphagenetwork/data/02_samtools_analysis/Beijing'
export GENOME_PATH=$genome_path
echo "The path to the genome data has been set to $GENOME_PATH."

# Check whether the index exists
if [ ! -f "./bj_sbatch.list" ]
then
    echo "The index does not exist, creating it ..."
    touch ./bj_sbatch.list
else
    echo "The index already exists, do you want to overwrite it? (y/n)"
    read answer
    if [ $answer == "y" ] || [ $answer == "Y" ]
    then
        echo "Overwriting the index ..."
        rm ./bj_sbatch.list
        touch ./bj_sbatch.list
    else
        echo "The index will not be overwritten."
        exit
    fi
fi

# Add the genome data to the config file "samples.index"
echo "Adding the genome data to the config file..."
for file_fq1 in $(ls ${GENOME_PATH})
do
    echo "$GENOME_PATH/$file_fq1" >> ./bj_sbatch.list
    sample_name=$(echo "$file_fq1" | grep -oE 'BJ[0-9]{3}')
    echo "The genome data of $sample_name has been added to the config file."
done

# Check if any error occurred during the process
if [ $? -eq 0 ]; then
    echo "Success!"
else
    echo "Warning: An error occurred during the process."
    exit 1
fi
