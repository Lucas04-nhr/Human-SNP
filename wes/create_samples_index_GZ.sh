#! /bin/bash

# Set the path to the genome data
echo "Please enter the path to the genome data of Guangzhou:"
read genome_path
# The path of the genome data is '/mnt/raid6/bacphagenetwork/data/samtools_results/Guangzhou'
# First check whether the directory exists.
if [ ! -d $genome_path ]
then
    echo "The directory does not exist, please check the path."
    exit
fi

# Set the path to the genome data
export GENOME_PATH=$genome_path
echo "The path to the genome data has been set to $GENOME_PATH."

# Check whether the index exists
if [ ! -f "./Guangzhou/samples.index" ]
then
    echo "The index does not exist, creating it ..."
    if [ -d ./Guangzhou ]
    then
        echo "The directory already exists."
    else
        echo "The directory does not exist, creating it ..."
        mkdir ./Guangzhou
    fi

else
    echo -c "The index of Guangzhou already exists, do you want to overwrite it? (y/n) "
    read answer
    if [ $answer == "y" ] || [ $answer == "Y" ]
    then
        echo "Overwriting the index ..."
        rm -rf ./Guangzhou/samples.index
    else
        echo "The index will not be overwritten."
        exit
    fi
fi

# Add the genome data to the config file "samples.index"
echo "Adding the genome data to the config file..."
for file in $GENOME_PATH/*.bam
do
    file_name=$(basename $file)
    sample_name=$(echo "$file_name" | grep -oE 'BJ[0-9]{3}')
    echo "$file_name    $file   0" >> ./Guangzhou/samples.index
done



# Check if any error occurred during the process
if [ $? -eq 0 ]; then
    echo "Success!"
else
    echo "Warning: An error occurred during the process."
    exit 1
fi