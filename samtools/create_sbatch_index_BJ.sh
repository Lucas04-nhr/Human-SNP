#! /bin/bash

# Set the path to the genome data
echo "Please enter the path to the processed *.sam data:"
read genome_path
# The path of the processed *.sam data is '/mnt/raid6/bacphagenetwork/data/bwa_analysis/Beijing'
export GENOME_PATH=$genome_path
echo "The path to the genome data has been set to $GENOME_PATH."

# Check whether the index exists
if [ ! -f "./bj_sbatch.list" ]
then
    echo "The index does not exist, creating it ..."

else
    echo "The index already exists, do you want to overwrite it? (y/n)"
    read answer
    if [ $answer == "y" ] || [ $answer == "Y" ]
    then
        echo "Overwriting the index ..."
    else
        echo "The index will not be overwritten."
        exit
    fi
fi

# Add the genome data to the config file "samples.index"
echo "Adding the genome data to the config file..."
for file_fq1 in $(ls ${GENOME_PATH})
do
    echo "$file_fq1" >> ./slurm/bj_sbatch.list
    echo "The genome data of $sample_name has been added to the config file."
done

# Check if any error occurred during the process
if [ $? -eq 0 ]; then
    echo "Success!"
else
    echo "Warning: An error occurred during the process."
    exit 1
fi
