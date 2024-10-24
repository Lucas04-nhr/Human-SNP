#! /bin/bash

# Set the path to the genome data
echo "Please enter the path to the *.g.vcf.gz files of Beijing:"
read genome_path
# The path of the genome data is '/mnt/raid6/bacphagenetwork/data/07_HaplotypeCaller/Beijing'
export GENOME_PATH=$genome_path
echo "The path to the genome data has been set to $GENOME_PATH."

# Print the test file range
echo "THIS SCRIPT IS FOR TESTING PURPOSES ONLY!"
echo "PRESS CTRL+C TO ABORT THE PROCESS!"
sleep 5
echo "The test file range is BJ001~BJ009."

# Check whether the index exists
if [ ! -f "./BJ_sbatch.list" ]
then
    echo "The index does not exist, creating it ..."

else
    echo "The index already exists, do you want to overwrite it? (y/n)"
    read answer
    if [ $answer == "y" ] || [ $answer == "Y" ]
    then
        echo "Overwriting the index ..."
        rm BJ_sbatch.list
    else
        echo "The index will not be overwritten."
        exit
    fi
fi

# Add the genome data to the config file "samples.index"
echo "Adding the genome data to the config file..."
for file in $GENOME_PATH/*.vcf.gz
do
    sample_name=$(basename "$file" .vcf.gz)
    if [[ $sample_name == BJ00* ]]; then
        echo "Sample name ${sample_name} is in the test file range."
        echo "Processing ${sample_name}..."
        echo $file >> BJ_sbatch.list
    else
        echo "Sample name ${sample_name} is not in the test file range."
        continue
    fi
done

# Check if any error occurred during the process
if [ $? -eq 0 ]; then
    echo "Success!"
else
    echo "Warning: An error occurred during the process."
    exit 1
fi