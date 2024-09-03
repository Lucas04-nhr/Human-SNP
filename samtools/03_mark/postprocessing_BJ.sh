#! /bin/bash

# Initialize the environment
export INPUT_PATH='/mnt/raid6/bacphagenetwork/data/03_samtools_sorted/Beijing'
export OUTPUT_PATH='/mnt/raid6/bacphagenetwork/data/04_dulplicate_marked/Beijing'

# Copy the metrics file
echo "Copying the metrics file..."

# Check the directory
if [ ! -d "./metrics/Beijing" ]
then
    echo "Creating the directory..."
    mkdir ./metrics
    mkdir ./metrics/Beijing
fi

for file in $(ls $OUTPUT_PATH/*.marked.metrics.txt)
do
    echo "Copying $file..."
    cp $file ./metrics/Beijing \
    || { echo "Error: Failed to copy $file"; exit 1; }
done

echo "The metrics files of Beijing have been copied."