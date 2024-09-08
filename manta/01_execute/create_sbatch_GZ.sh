#! /bin/bash

# Check whether the index exists
if [ ! -f "./GZ_sbatch.list" ]
then
    echo "The index does not exist, creating it ..."
    touch ./GZ_sbatch.list
else
    echo "The index already exists, do you want to overwrite it? (y/n)"
    read answer
    if [ $answer == "y" ] || [ $answer == "Y" ]
    then
        echo "Overwriting the index ..."
        rm ./GZ_sbatch.list
        touch ./GZ_sbatch.list
    else
        echo "The index will not be overwritten."
        exit
    fi
fi

# Set the path to the genome data
echo "Please enter the path to the parent execute folder of Guangzhou:"
read genome_path
# The path of the original *.bam data is '/mnt/raid6/bacphagenetwork/data/07_manta/01_exec/Guangzhou'
export GENOME_PATH=$genome_path
echo "The path to the genome data has been set to $GENOME_PATH."

# Add the genome data to the config file "samples.index"
echo "Adding the execute folders to the config file..."
for folder in $(ls ${GENOME_PATH})
do
    echo "$GENOME_PATH/$folder" >> ./GZ_sbatch.list
    sample_name=$(echo "$folder" | grep -oE 'GZ[0-9]{3}')
    echo "The execute folder of $sample_name has been added to the config file."
done

# Create the log folder
if [ ! -d "./log/Guangzhou" ]
then
    mkdir -p ./log/Guangzhou
    echo "The log folder has been created."
else
    echo "The log folder already exists."
fi

# Check if any error occurred during the process
if [ $? -eq 0 ]; then
    echo "Success!"
else
    echo "Warning: An error occurred during the process."
    exit 1
fi
