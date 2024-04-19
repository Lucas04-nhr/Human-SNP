#! /bin/bash

# Check whether conda is installed
if ! command -v conda &> /dev/null
then
    echo "Conda could not be found, please install it first."
    exit
fi

# Check whether the environment exists
if conda env list | grep -q "wescall"
then
    echo "Great! The environment already exists."
    # Activate the environment
    echo "Activating the environment..."
    conda activate wescall
else
    echo "Creating the environment..."
    conda env create -n wescall -f requirements.txt
    echo "The environment has been created, activating it..."
    conda activate wescall
fi

# Create the samples index
bash create_samples_index.sh
if [ $? -ne 0 ]; then
    echo "Failed to create the samples index. Please check the error message."
fi

# Process the samples
bash processing.sh
if [ $? -ne 0 ]; then
    echo "Failed to process the samples. Please check the error message."
fi