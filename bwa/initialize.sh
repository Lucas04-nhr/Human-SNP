#! /bin/bash

# Check whether conda is installed
if ! command -v conda &> /dev/null
then
    echo "Conda could not be found, please install it first."
    exit
fi

# Check whether the environment exists
if conda env list | grep -q "bwa"
then
    echo "Great! The environment already exists."
else
    echo "Creating the environment..."
    conda env create -f requirements.txt
    echo "The environment has been created."
fi

