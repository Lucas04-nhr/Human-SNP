#! /bin/bash

# Initialize the environment
echo "Initializing..."
export INPUT_PATH='/mnt/raid6/bacphagenetwork/data/03_samtools_sorted/Beijing'
export OUTPUT_PATH='/mnt/raid6/bacphagenetwork/data/04_dulplicate_marked/Beijing'

# Check the metrics file
echo "Checking the metrics file..."
if [ ! -f "$OUTPUT_PATH/*.marked.metrics.txt" ]
then
    echo "Error: $OUTPUT_PATH/*.marked.metrics.txt does not exist."
    exit 1
else
    echo "The file $OUTPUT_PATH/*.marked.metrics.txt exists."
fi

echo "Initializing complete."

# Copy the metrics file
echo "Copying the metrics file..."

# Check the directory
if [ ! -d "./metrics/Beijing" ]
then
    echo "Creating the directory..."
    mkdir ./metrics
    mkdir ./metrics/Beijing
fi

cp $OUTPUT_PATH/*.marked.metrics.txt ./metrics/Beijing || { echo "Error: Failed to copy the metrics file of $OUTPUT_PATH/*.marked.metrics.txt"; exit 1; }

echo "The metrics file has been copied."