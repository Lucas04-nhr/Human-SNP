#! /bin/bash

# Set the path to the genome data
echo "Please enter the path to the sorted *bam files of Guangzhou:"
read genome_path
# The path of the genome data is '/mnt/raid6/bacphagenetwork/data/03_sort/Guangzhou'
export GENOME_PATH=$genome_path
echo "The path to the genome data has been set to $GENOME_PATH."

# Check whether the index exists
if [ ! -f "./GZ_sbatch.list" ]
then
    echo "The index does not exist, creating it ..."

else
    echo "The index already exists, do you want to overwrite it? (y/n)"
    read answer
    if [ $answer == "y" ] || [ $answer == "Y" ]
    then
        echo "Overwriting the index ..."
        rm GZ_sbatch.list
    else
        echo "The index will not be overwritten."
        exit
    fi
fi

# Add the genome data to the config file "samples.index"
echo "Adding the genome data to the config file..."
for file in $GENOME_PATH/*.bam
do
    sample_name=$(basename "$file" .sorted.bam)
    echo "Processing ${sample_name}..."
    echo $file >> GZ_sbatch.list
done

# Check if any error occurred during the process
if [ $? -eq 0 ]; then
    echo "Success!"
else
    echo "Warning: An error occurred during the process."
    exit 1
fi