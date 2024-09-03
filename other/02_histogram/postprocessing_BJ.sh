#! /bin/bash

# Initialize the environment
export OUTPUT_PATH='/mnt/raid6/bacphagenetwork/data/07_histogram/Beijing'

# Copy the metrics file
echo "Copying the metrics file..."

# Check the directory
if [ ! -d "./metrics/Beijing" ]
then
    echo "Creating the directory..."
    if [ ! -d "./metrics" ]
    then
        mkdir ./metrics
    fi
    mkdir ./metrics/Beijing
    mkdir ./metrics/Beijing/pdf
fi

for file in $(ls $OUTPUT_PATH/*.histogram.txt)
do
    export sample_name=$(echo "$file" | grep -oE 'BJ[0-9]{3}')
    echo "Copying metrics file of $sample_name..."
    cp $file ./metrics/Beijing \
    || { echo "Error: Failed to copy $file"; exit 1; }
done

for file in $(ls $OUTPUT_PATH/*.histogram.pdf)
do
    export sample_name=$(echo "$file" | grep -oE 'BJ[0-9]{3}')
    echo "Copying histogram of $sample_name..."
    cp $file ./metrics/Beijing/pdf \
    || { echo "Error: Failed to copy $file"; exit 1; }
done

echo "The metrics files of Beijing have been copied."