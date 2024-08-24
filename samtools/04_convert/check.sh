#! /bin/bash

# Check if the input folder exists

echo "Checking the input folder..."

if [ ! -d "/mnt/raid6/bacphagenetwork/data/05_format_converted" ]
then
    echo "The input folder does not exist, please check manually ..."
    exit 1
else
    echo "The input folder already exists."
fi

# Create the log folder

echo "Creating the log folder..."

if [ -d "./log" ]
then
    echo "The log folder already exists, do you want to back up the log folder? (y/n)"
    read answer
    if [ $answer == "y" ] || [ $answer == "Y" ]
    then
        echo "Backing up the log folder ..."
        mv log log_bak_$(date +%Y%m%d%H%M%S)
    else
        echo "The original log folder will be removed."
        rm -rf log
    fi
fi

mkdir log
mkdir log/Beijing
mkdir log/Guangzhou
echo "The log folder has been created."

# Create the output folder

echo "Creating the output folder..."

if [ -d "/mnt/raid6/bacphagenetwork/data/05_format_converted" ]
then
    echo "The output folder already exists, do you want to back up the output folder? (y/n)"
    read answer
    if [ $answer == "y" ] || [ $answer == "Y" ]
    then
        echo "Backing up the output folder ..."
        mv /mnt/raid6/bacphagenetwork/data/05_format_converted /mnt/raid6/bacphagenetwork/data/05_format_converted_bak_$(date +%Y%m%d%H%M%S)
    else
        echo "The original output folder will be removed."
        rm -rf /mnt/raid6/bacphagenetwork/data/05_format_converted
    fi
fi

mkdir /mnt/raid6/bacphagenetwork/data/05_format_converted
mkdir /mnt/raid6/bacphagenetwork/data/05_format_converted/Beijing
mkdir /mnt/raid6/bacphagenetwork/data/05_format_converted/Guangzhou

echo "The output folder has been created."

# Ask for creating sbatch files

echo "Do you want to create sbatch files? (y/n)"
read answer
if [ $answer != "y" ] && [ $answer != "Y" ]
then
    echo "The process has been terminated."
    exit 0
fi

# Creating sbatch files

echo "Creating sbatch files ..."

bash ./create_sbatch_BJ.sh
bash ./create_sbatch_GZ.sh

echo "The sbatch files have been created."