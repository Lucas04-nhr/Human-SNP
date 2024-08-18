#! /bin/bash

# Initialize the environment
export INPUT_PATH='/mnt/raid6/bacphagenetwork/data/03_samtools_sorted/Guangzhou'
export OUTPUT_PATH='/mnt/raid6/bacphagenetwork/data/04_dulplicate_marked/Guangzhou'

# Copy the metrics file
echo "Copying the metrics file..."

# Check the directory
if [ ! -d "./metrics" ]
then
    mkdir ./metrics
fi

if [ ! -d "./metrics/Guangzhou" ]
then
    echo "Creating the directory..."
    mkdir ./metrics/Guangzhou
fi

for file in $(ls $OUTPUT_PATH/*.marked.metrics.txt)
do
    echo "Copying $file..."
    cp $file ./metrics/Guangzhou \
    || { echo "Error: Failed to copy $file"; exit 1; }
done

echo "The metrics files of Guangzhou have been copied."